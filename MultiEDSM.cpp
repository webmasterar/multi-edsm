/*
    MultiEDSM: Multiple Elastic Degenerate String Matching

    Copyright (C) 2017 Ahmad Retha and Solon P. Pissis.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <cstdlib>
#include <climits>
#include <string>
#include <vector>
#include <queue>
#include <ctime>
#include <map>
#include <divsufsort64.h>
#include <sdsl/util.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_trees.hpp>
#include "MyUMSA.hpp"
#include "MultiEDSM.hpp"

using namespace sdsl;
using namespace std;

/**
 * @constructor
 * @param alphabet Alphabet used by the patterns
 * @param patterns The patterns to look for
 * @param maxNoBitVectorsStorable The maximum number of patterns to store in memory as uncompressed bitvectors
 */
MultiEDSM::MultiEDSM(const string & alphabet, const vector<string> & patterns, const unsigned int maxNoBitVectorsStorable)
{
    this->alphabet = alphabet;
    if (patterns.size() == 0) {
        cerr << "Error: Empty pattern set!" << endl;
        throw 1;
    }
    this->maxK = maxNoBitVectorsStorable;
    this->minP = UINT_MAX;
    this->maxP = 0;
    this->f = 0;
    this->F = 0;
    this->d = 0;
    this->D = 0;
    this->M = 0;
    this->R = 0;
    this->Np = 0;
    this->Nm = 0;
    this->pos = 0;
    this->primed = false;
    this->duration = 0;
    this->reportOnce = false;
    this->reportPatterns = true;
    this->preprocessPatterns(patterns);
}

/**
 * @destructor
 */
MultiEDSM::~MultiEDSM()
{
    delete this->umsa;
}

/**
 * Proprocess the patterns - initialize the bitvector and OccVector tools.
 * Complexity:
 *  - Shift-And time and space O(M + sigma*[M/w]).
 *  - Suffix Tree time and space O(M + k), where k is number of patterns.
 *  - OccVector construction time O(M + k), space O((M + k) * [M/w])
 *
 * @param patterns
 */
void MultiEDSM::preprocessPatterns(const vector<string> & patterns)
{
    //start timer
    clock_t start = clock();

    //initialize Shift-And searching tool
    cout << "UMSA init..." << endl;
    this->umsa = new MyUMSA(this->alphabet, this->reportPatterns);

    //add patterns to bitvector and build suffix tree string
    string p = "";
    unsigned int h, i, j;
    for (i = 0; i < patterns.size(); i++)
    {
        h = patterns[i].length();
        if (h == 0) {
            cerr << "Error: zero-length pattern given!" << endl;
            throw 1;
        }

        //create bitvector of pattern - for Shift-And
        this->umsa->addPattern(patterns[i]);

        //add pattern to p and append with seperator (#) for building suffix tree of patterns
        p += patterns[i];
        p += SEPARATOR_CHAR;

        //update STpIdx2BVIdx datastructures with correct index of STp match for
        //bitvector and pattern id based on bit position in bitvector
        this->M += h;
        for (j = 0; j < h; j++)
        {
            if (this->alphabet.find(patterns[i][j]) == string::npos)
            {
                cerr << "Error: Invalid character found in pattern: " << (int) patterns[i][j] << endl;
                throw 1;
            }
            this->STpIdx2BVIdx.push_back(this->R + j - i);
        }
        this->STpIdx2BVIdx.push_back(SEPARATOR_DIGIT);
        this->R += h + 1;

        //find min/max pattern length
        if (h < this->minP) {
            this->minP = h;
        }
        if (h > this->maxP) {
            this->maxP = h;
        }
    }
    p.pop_back(); //remove unnecessary separator char as \0 already tagged on end of c_str

    //tool to get pattern id from any given position in the bitvector/p
    cout << "Pos2Pat..." << endl;
    this->Pos2PatId = this->umsa->getPatternPositions();

    //construct the suffix tree of patterns with seperators
    cout << "SuffixTree..." << endl;
    construct_im(this->STp, p.c_str(), sizeof(char));

    //construct occVector datastructure
    cout << "OccVector..." << endl;
    this->constructOV7(p);

    //stop timer
    this->duration += clock() - start;

    cout << "Preproccessing " << patterns.size() << " patterns of size M=" \
         << this->M << " completed in " << this->getDuration() << "s." << endl << endl;
}

/**
 * - Minimalist + pruning.
 * - Uses word vector array OVMem.
 * - Recursive tree traversal.
 *
 * Construct the occVector tool. Requires O(M + k) time (i.e. the number of nodes
 * in a suffix tree is at most 2M, and k is the number of patterns, represented
 * by different seperators). Requires O((M + k) * [(M + k) / w]) space, because
 * each node stores at most O([(M + k) / w]) computer words.
 */
void MultiEDSM::constructOV()
{
    unsigned int i, numNodes = this->STp.nodes();
    cout << "NumNodes: " << numNodes << endl;
    for (i = 0; i < numNodes; i++) {
        WordVector v(1, 0ul);
        this->OVMem.push_back(v);
    }
    cout << "RecAssign" << endl;
    this->recAssignOVMem(this->STp.root());
}
/**
 * @deprecated This method uses O((M+k)/w * (M+k)) space which is a lot more than
 * MultiEDSM::recAssignOVMemRestricted() which replaces it.
 *
 * Recursively traverses the suffix tree STp and denotes starting positions of
 * substrings, encoding them into WordVectors, and combining them as it moves
 * back up the tree. Space efficient because it only uses as many computer words
 * as is necessary.
 */
WordVector MultiEDSM::recAssignOVMem(const cst_node_t & u)
{
    unsigned int id = this->STp.id(u);

    if (this->STp.is_leaf(u))
    {
        unsigned int h = this->STp.sn(u); //h is index in suffix tree string

        //make sure it is not first letter or last letter in a pattern
        if (h > 0 && this->STpIdx2BVIdx[h - 1] != SEPARATOR_DIGIT && \
           (h + 1) < ((int)this->R - 1) && this->STpIdx2BVIdx[h + 1] != SEPARATOR_DIGIT)
        {
            //calculate bit position
            long long int i = this->STpIdx2BVIdx[h] - 1;    //i is index in bitvector
            unsigned int j = this->OVMem[id].size();        //current size of the OVMem[id] WordVector
            unsigned int wordIdx = (unsigned int) ((double)i / (double)BITSINWORD);
            unsigned int bitShifts = i % BITSINWORD;

            //ensure word vector has enough words in it
            while (j <= wordIdx) {
                this->OVMem[id].push_back(0ul);
                j++;
            }

            //write position to OVMem[id]
            this->OVMem[id][wordIdx] = this->OVMem[id][wordIdx] | (1ul << bitShifts);
        }
    }
    else
    {
        for (const auto & child : this->STp.children(u)) {
            this->WordVectorOR_IP(this->OVMem[id], this->recAssignOVMem(child));
        }
    }

    return this->OVMem[id];
}

/**
 * - Recursive tree traversal with maps.
 * - Uses map for node ids.
 * - Recursive tree traversal with depth checking.
 *
 * Construct the occVector tool. Requires O(M + k) time (i.e. the number of nodes
 * in a suffix tree is at most 2M, and k is the number of patterns, represented
 * by different seperators). Requires O((M + k) * [(M + k) / w]) space, because
 * each node stores at most O([(M + k) / w]) computer words.
 */
void MultiEDSM::constructOV2()
{
    unsigned int i, numNodes = this->STp.nodes();
    cout << "NumNodes: " << numNodes << endl;
    cout << "RecAssign" << endl;
    this->recAssignOVMem2(this->STp.root(), 0);
}
WordVector MultiEDSM::recAssignOVMem2(const cst_node_t & u, const unsigned int currDepth)
{
    unsigned int id = this->STp.id(u);
    unsigned int depth = (currDepth <= (this->maxP - 2)) ? this->STp.node_depth(u) : currDepth;

    if (this->STp.is_leaf(u))
    {
        unsigned int h = this->STp.sn(u); //h is index in suffix tree string

        //make sure it is not first letter or last letter in a pattern
        if (
            (h > 0 && ((h + 1) < (this->R - 1))) && \
            (this->STpIdx2BVIdx[h-1] != SEPARATOR_DIGIT) && \
            (this->STpIdx2BVIdx[h]   != SEPARATOR_DIGIT) && \
            (this->STpIdx2BVIdx[h+1] != SEPARATOR_DIGIT)
        )
        {
            //calculate bit position
            unsigned int i = this->STpIdx2BVIdx[h] - 1; //i is index in bitvector
            unsigned int wordIdx = (unsigned int) ((double)i / (double)BITSINWORD);
            unsigned int bitShifts = i % BITSINWORD;

            if (depth > (this->maxP - 2))
            {
                WordVector v(wordIdx + 1, 0ul);
                v[wordIdx] = v[wordIdx] | (1ul << bitShifts);
                return v;
            }
            else
            {
                if (this->OVMem2.find(id) == this->OVMem2.end())
                {
                    WordVector v(wordIdx + 1, 0ul);
                    this->OVMem2[id] = v;
                }

                //write position to OVMem2[id]
                this->OVMem2[id][wordIdx] = this->OVMem2[id][wordIdx] | (1ul << bitShifts);
                return this->OVMem2[id];
            }
        }
        else
        {
            WordVector v;
            return v;
        }
    }
    else
    {
        //now combine children with self
        if (depth <= (this->maxP - 2))
        {
            if (this->OVMem2.find(id) == this->OVMem2.end())
            {
                WordVector v;
                this->OVMem2[id] = v;
            }
            for (const auto & child : this->STp.children(u)) {
                this->WordVectorOR_IP(this->OVMem2[id], this->recAssignOVMem2(child, depth + 1));
            }
            return this->OVMem2[id];
        }
        else
        {
            WordVector v;
            for (const auto & child : this->STp.children(u)) {
                this->WordVectorOR_IP(v, this->recAssignOVMem2(child, depth + 1));
            }
            return v;
        }
    }
}

/**
 * - Reverse level order tree traversal.
 * - Uses map for nodes.
 */
void MultiEDSM::constructOV3()
{
    //first create array of the levels that each node is in in the suffix tree by
    //traversing the tree in level order. We need to do this because SDSL takes
    //O(depth) time for each node depth query but we want O(1) time. To get this
    //information we require time and extra space O(n + depth).
    cout << "Creating nodeLevels" << endl;
    vector<vector<unsigned int>> nodeLevels;
    queue<cst_node_t> q;
    q.push(this->STp.root());
    unsigned int level = 1;     //current level
    unsigned int dc = 1;        //decrement counter
    unsigned int ic = 0;        //increment counter
    cst_node_t currNode;
    while (!q.empty())
    {
        currNode = q.front();
        q.pop();
        for (const auto & child : this->STp.children(currNode))
        {
            if (nodeLevels.size() <= level) {
                vector<unsigned int> v;
                v.push_back(this->STp.id(child));
                nodeLevels.push_back(v);
            } else {
                nodeLevels[level].push_back(this->STp.id(child));
            }
            q.push(child);
            ic++;
        }
        dc--;
        if (dc == 0) {
            level++;
            dc = ic;
            ic = 0;
        }
    }
    cout << "Done creating nodeLevels" << endl;

    //Create OVMem then loop through the nodes in reverse level order, calculating
    //the bitvector and storing it in OVMem. Wipes memory as it rolls back up
    //towards the root
    int j;
    unsigned int nodeCreationCounter, leafCreationCounter, totalNodesCounter, totalLeafSize;
    unsigned int parent_id, root_id = this->STp.id(this->STp.root());
    unsigned int i, numNodes = this->STp.nodes();
    // cout << "Initializing OVMem" << endl;
    // sleep(5);
    // for (i = 0; i < numNodes; i++) {
    //     WordVector v(1, 0ul);
    //     v.shrink_to_fit();
    //     this->OVMem.push_back(v);
    // }
    // cout << "Done!" << endl;
    // sleep(5);
    cout << "Tree climbing started" << endl;
    for (j = nodeLevels.size() - 1; j >= 0; j--)
    {
        nodeCreationCounter = 0;
        leafCreationCounter = 0;
        totalNodesCounter = 0;
        totalLeafSize = 0;

        for (unsigned int id : nodeLevels[j])
        {
            currNode = this->STp.inv_id(id);

            if (id != root_id)
            {
                parent_id = this->STp.id(this->STp.parent(currNode));

                if (this->STp.is_leaf(currNode))
                {
                    unsigned int h = this->STp.sn(currNode);  //h is index in suffix tree string

                    //make sure it is not first letter or last letter in a pattern or seperator char
                    if (
                        (h > 0 && ((h + 1) < (this->R - 1))) && \
                        (this->STpIdx2BVIdx[h-1] != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[h]   != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[h+1] != SEPARATOR_DIGIT)
                    ) {
                        //calculate bit position
                        long long int k = this->STpIdx2BVIdx[h] - 1;  //i is index in bitvector
                        unsigned int l = 1;                           //current size of a
                        unsigned int wordIdx = (unsigned int) ((double)k / (double)BITSINWORD);
                        unsigned int bitShifts = k % BITSINWORD;

                        //create word vector with enough words in it
                        WordVector a(wordIdx + 1, 0ul);

                        //write position to a, OVMem[id] if applicable
                        a[wordIdx] = a[wordIdx] | (1ul << bitShifts);
                        if (j <= this->maxP) {
                            this->OVMem3[id] = a;
                        }

                        //create the parent with a if not exists or combine a with parent if it does
                        if (this->OVMem3.find(parent_id) == this->OVMem3.end()) {
                            this->OVMem3[parent_id] = a;
                            leafCreationCounter++;
                            totalLeafSize += wordIdx + 1;
                        } else {
                            this->WordVectorOR_IP(this->OVMem3[parent_id], a);
                        }
                    }
                }
                else
                {
                    //if node already exists
                    if (this->OVMem3.find(id) != this->OVMem3.end())
                    {
                        if (this->OVMem3.find(parent_id) == this->OVMem3.end()) {
                            this->OVMem3[parent_id] = this->OVMem3[id];
                            nodeCreationCounter++;
                        } else {
                            this->WordVectorOR_IP(this->OVMem3[parent_id], this->OVMem2[id]);
                        }
                        if (this->STp.depth(currNode) >= this->maxP) {
                            // this->OVMem[id].clear();
                            this->OVMem3.erase(id);
                        }
                    }
                }
            }
        }

        totalNodesCounter = nodeLevels[j].size();

        cout << "Level " << j << " contains " << totalNodesCounter << " nodes. " \
             << leafCreationCounter << " leaves and " << nodeCreationCounter \
             << " nodes were created. Total size: " << (totalLeafSize * sizeof(WORD))<< endl;

        nodeLevels[j].clear();
    }
}

/**
 * - Minimalist.
 * - Uses WordVector Array OVMem.
 * - Recursive tree traversal.
 */
void MultiEDSM::constructOV4()
{
    unsigned int i, numNodes = this->STp.nodes();
    cout << "NumNodes: " << numNodes << endl;
    for (i = 0; i < numNodes; i++) {
        WordVector v(1, 0ul);
        v.shrink_to_fit();
        this->OVMem4.push_back(v);
    }
    this->recAssignOVMem4(this->STp.root());
}
WordVector & MultiEDSM::recAssignOVMem4(const cst_node_t & u)
{
    unsigned int id = this->STp.id(u);

    if (this->STp.is_leaf(u))
    {
        unsigned int sn = this->STp.sn(u);
        long long int us = (long long int)this->STpIdx2BVIdx[sn] - 1;
        if (us > SEPARATOR_DIGIT)
        {
            unsigned int wordIdx = (unsigned int) ((double)us / (double)BITSINWORD);
            unsigned int bitShifts = (unsigned int) us % BITSINWORD;

            WordVector a(wordIdx + 1, 0ul);
            a.shrink_to_fit();
            a[wordIdx] = a[wordIdx] | (1ul << bitShifts);
            this->OVMem4[id] = a;
        }
    }
    else
    {
        for (const auto & child : this->STp.children(u)) {
            this->WordVectorOR_IP(this->OVMem4[id], this->recAssignOVMem4(child));
        }
    }

    return this->OVMem4[id];
}

/**
 * @deprecated idea abandoned
 * - Complex with SDSL types.
 * - Storage in two containers - OVMem5 map and OVMemLeaves map.
 * - Suffix Array traversal and level order tree traversal.
 */
void MultiEDSM::constructOV5()
{
    //
    // Step 1: Fill in OVMemLeaves -- abandoned due to memory problems
    //

    //store leaves in OVMemLeaves
    // cout << "Handle leaves..." << endl;
    // long long int k, i = this->STp.csa.size() - 1; //last position
    // unsigned int parent_id, wordIdx, bitShifts;
    // while (i >= 0)
    // {
    //     k = this->STpIdx2BVIdx[this->STp.csa[i]];
    //     if (k != SEPARATOR_DIGIT)
    //     {
            //make leaf
            // bit_vector bv(k, 0, 1);
            // k = k - 1;
            // bv[k] = 1;
            // RRR r(bv);
            // this->OVMemLeaves[i] = r;

            //then store WordVector in parent -- this is problematic as there are as many leaves as the size of the input O(M) * O(M/w) = O(M^2/w)
            // wordIdx = (unsigned int) ((double)k / (double)BITSINWORD);
            // bitShifts = (unsigned) (k % BITSINWORD);
            // WordVector wv(wordIdx + 1, 0ul);
            // wv[wordIdx] = wv[wordIdx] | (1ul << bitShifts);
            // parent_id = this->STp.id(this->STp.parent(this->STp.inv_id(i)));
            // if (this->OVMem5.find(parent_id) == this->OVMem5.end()) {
            //     this->OVMem5[parent_id] = wv;
            // } else {
            //     this->WordVectorOR_IP(this->OVMem5[parent_id], wv);
            // }
    //     }
    //     i--;
    // }
}

void MultiEDSM::constructOV6(const string & p)
{
    //
    // Step 1: Traverse the tree level order down to maxP-2 (inclusive) levels
    // making sure not to follow edges beginning with a separator char
    //
    cout << "Creating nodeLevels..." << endl;
    vector<vector<unsigned int>> nodeLevels;
    vector<unsigned int> v;
    v.push_back(this->STp.id(this->STp.root()));
    nodeLevels.push_back(v);
    queue<cst_node_t> q;
    q.push(this->STp.root());
    unsigned int id, maxId = 0;
    unsigned int maxDepth = (unsigned int) max((int)1, (int)this->maxP - 2); //max level to traverse in level order
    unsigned int level = 1;     //current level -- root is level 0, but we are inserting children into q, so level = currNodeLevel+1
    unsigned int dc = 1;        //decrement counter
    unsigned int ic = 0;        //increment counter
    cst_node_t currNode;
    while (!q.empty() && level <= maxDepth)
    {
        currNode = q.front();
        q.pop();
        for (const auto & child : this->STp.children(currNode))
        {
            // id = this->STp.id(child);
            // if (nodeLevels.size() <= level) {
            //     vector<unsigned int> v;
            //     v.push_back(id);
            //     nodeLevels.push_back(v);
            // } else {
            //     nodeLevels[level].push_back(id);
            // }
            // q.push(child);
            // ic++;

            //or

            if (!this->STp.is_leaf(child))
            {
                unsigned int lb = this->STp.lb(child);
                unsigned int sn = this->STp.csa[lb];
                char nextChar = p[sn + level - 1];
                if (nextChar != SEPARATOR_CHAR) {
                    id = this->STp.id(child);
                    if (id > maxId) {
                        maxId = id;
                    }
                    if (nodeLevels.size() <= level) {
                        vector<unsigned int> v;
                        v.push_back(id);
                        nodeLevels.push_back(v);
                    } else {
                        nodeLevels[level].push_back(id);
                    }
                    q.push(child);
                    ic++;
                }
            }
            else
            {
                id = this->STp.id(child);
                if (nodeLevels.size() <= level) {
                    vector<unsigned int> v;
                    v.push_back(id);
                    nodeLevels.push_back(v);
                } else {
                    if (id > maxId) {
                        maxId = id;
                    }
                    nodeLevels[level].push_back(id);
                }
                q.push(child);
                ic++;
            }
        }
        dc--;
        if (dc == 0) {
            level++;
            dc = ic;
            ic = 0;
        }
    }

    for (unsigned int x = 0; x < nodeLevels.size(); x++)
    {
        cout << "Level " << x << " contains " << nodeLevels[x].size() << " nodes." << endl;
    }

    // With separator char checking:
    // Level 0 contains 1 nodes.
    // Level 1 contains 5 nodes.
    // Level 2 contains 17 nodes.
    // Level 3 contains 65 nodes.
    // Level 4 contains 257 nodes.
    // Level 5 contains 1025 nodes.
    // Level 6 contains 4097 nodes.
    // Level 7 contains 16385 nodes.
    // Level 8 contains 65537 nodes.
    // Level 9 contains 260952 nodes.
    // Without separator char checking:
    // Level 0 contains 1 nodes.
    // Level 1 contains 6 nodes.
    // Level 2 contains 25 nodes.
    // Level 3 contains 113 nodes.
    // Level 4 contains 513 nodes.
    // Level 5 contains 2305 nodes.
    // Level 6 contains 10241 nodes.
    // Level 7 contains 45057 nodes.
    // Level 8 contains 196609 nodes.
    // Level 9 contains 849631 nodes.


    //
    // Step 2: For each node at maxDepth, recursively traverse its children and
    // encode the leaves to it in OVMem.
    //
    cout << "Initializing OVMem6 with " << (maxId + 1) << " nodes..." << endl;
    unsigned int i, numNodes = maxId + 1; //this->STp.nodes();
    this->OVMem6.reserve(numNodes);
    for (i = 0; i < numNodes; i++)
    {
        WordVector v; //(1, 0ul);
        // v.shrink_to_fit();
        this->OVMem6.push_back(v);
    }
    for (unsigned int nodeId : nodeLevels[maxDepth])
    {
        currNode = this->STp.inv_id(nodeId);
        if (!this->STp.is_leaf(currNode)) {
            cout << "Leaf searching " << nodeId;
            // this->recEncodeLeavesToNode(currNode, nodeId);
            // this->nonrecEncodeLeavesToNode(currNode, nodeId);
            this->csaEncodeLeavesToNode(currNode, nodeId);
            this->OVMem6[nodeId].shrink_to_fit();
            cout << ", size " << (unsigned int) ((double)(this->OVMem6[nodeId].size() * 8) / 1024.0 / 1024.0) << "MB" << endl;
        }
    }

    //
    // Step 3: Traverse the tree in reverse level order to encode children into
    // their parents from maxDepth down to level 1.
    //
    cout << "Tree climbing..." << endl;
    unsigned int parentId, h;
    for (i = maxDepth; i > 0; i--)
    {
        for (unsigned int nodeId : nodeLevels[i])
        {
            currNode = this->STp.inv_id(nodeId);
            parentId = this->STp.id(this->STp.parent(currNode));

            if (this->STp.is_leaf(currNode))
            {
                h = this->STp.sn(currNode);   //h is index in suffix tree string
                //make sure it is not first letter or last letter in a pattern or seperator char
                if (
                    (h > 0 && ((h + 1) < (this->R - 1))) && \
                    (this->STpIdx2BVIdx[h-1] != SEPARATOR_DIGIT) && \
                    (this->STpIdx2BVIdx[h]   != SEPARATOR_DIGIT) && \
                    (this->STpIdx2BVIdx[h+1] != SEPARATOR_DIGIT)
                ) {
                    //calculate bit position
                    long long int k = this->STpIdx2BVIdx[h] - 1;  //i is index in bitvector
                    unsigned int wordIdx = (unsigned int) ((double)k / (double)BITSINWORD);
                    unsigned int bitShifts = k % BITSINWORD;

                    //create word vector with enough words in it and assign to target node
                    WordVector a(wordIdx + 1, 0ul);
                    a[wordIdx] = a[wordIdx] | (1ul << bitShifts);
                    this->OVMem6[nodeId] = a;

                    //combine with parent
                    this->WordVectorOR_IP(this->OVMem6[parentId], this->OVMem6[nodeId]);
                }
            }
            else
            {
                this->WordVectorOR_IP(this->OVMem6[parentId], this->OVMem6[nodeId]);
            }
        }
    }
}
void MultiEDSM::recEncodeLeavesToNode(const cst_node_t & node, const unsigned int targetNodeId)
{
    if (this->STp.is_leaf(node))
    {
        unsigned int h = this->STp.sn(node);  //h is index in suffix tree string

        //make sure it is not first letter or last letter in a pattern or seperator char
        if (
            (h > 0 && ((h + 1) < (this->R - 1))) && \
            (this->STpIdx2BVIdx[h-1] != SEPARATOR_DIGIT) && \
            (this->STpIdx2BVIdx[h]   != SEPARATOR_DIGIT) && \
            (this->STpIdx2BVIdx[h+1] != SEPARATOR_DIGIT)
        ) {
            //calculate bit position
            long long int k = this->STpIdx2BVIdx[h] - 1;  //i is index in bitvector

            // unsigned int wordIdx = (unsigned int) ((double)k / (double)BITSINWORD);
            // unsigned int bitShifts = k % BITSINWORD;
            //
            // //create word vector with enough words in it and combine it with target node
            // WordVector a(wordIdx + 1, 0ul);
            // a[wordIdx] = a[wordIdx] | (1ul << bitShifts);
            // this->WordVectorOR_IP(this->OVMem6[targetNodeId], a);

            this->WordVectorSet1At(this->OVMem6[targetNodeId], (unsigned int)k);
        }
    }
    else
    {
        for (const auto & child : this->STp.children(node)) {
            this->recEncodeLeavesToNode(child, targetNodeId);
        }
    }
}
void MultiEDSM::nonrecEncodeLeavesToNode(const cst_node_t & node, const unsigned int targetNodeId)
{
    unsigned int h, k;
    queue<cst_node_t> q;
    q.push(node);
    cst_node_t currNode;
    while (!q.empty())
    {
        currNode = q.front();
        q.pop();
        for (const auto & child : this->STp.children(currNode))
        {
            if (this->STp.is_leaf(child))
            {
                h = this->STp.sn(child);  //h is index in suffix tree string

                //make sure it is not first letter or last letter in a pattern or seperator char
                if (
                    (h > 0 && ((h + 1) < (this->R - 1))) && \
                    (this->STpIdx2BVIdx[h-1] != SEPARATOR_DIGIT) && \
                    (this->STpIdx2BVIdx[h]   != SEPARATOR_DIGIT) && \
                    (this->STpIdx2BVIdx[h+1] != SEPARATOR_DIGIT)
                ) {
                    //calculate bit position
                    k = this->STpIdx2BVIdx[h] - 1;  //i is index in bitvector
                    this->WordVectorSet1At(this->OVMem6[targetNodeId], k);
                }
            }
            else
            {
                q.push(child);
            }
        }
    }
}
void MultiEDSM::csaEncodeLeavesToNode(const cst_node_t & node, const unsigned int targetNodeId)
{
    unsigned int lb = this->STp.lb(node);
    unsigned int rb = this->STp.rb(node);
    unsigned int h, i, k;
    for (i = lb; i < rb + 1; i++)
    {
        h = this->STp.csa[i];
        if (
            (h > 0 && ((h + 1) < (this->R - 1))) && \
            (this->STpIdx2BVIdx[h-1] != SEPARATOR_DIGIT) && \
            (this->STpIdx2BVIdx[h]   != SEPARATOR_DIGIT) && \
            (this->STpIdx2BVIdx[h+1] != SEPARATOR_DIGIT)
        ) {
            //calculate bit position
            k = this->STpIdx2BVIdx[h] - 1;  //i is index in bitvector
            this->WordVectorSet1At(this->OVMem6[targetNodeId], k);
        }
    }
    //
    // The plan is to use getMaxRangeInSet() and encodeCompressedVector() to
    // create compressed bitvectors in OVMem6 instead of uncompressed ones.
    // Uncompressed takes up [M/w] space which is 5MB for my 1/8 data sample,
    // requiring total of 1.5TB memory! Compressed bitvectors will take up far
    // less space.
    //
}
/**
 * Takes a set of integers and returns the largest difference or range between
 * the numbers. 0 to min(s) is also considered. Uses the linear time and space
 * Maximal Gap algorithm to achieve this. e.g. s = [9, 19, 27] --> 10 (19 - 9),
 * or s = [50, 60, 90] --> 50 (50 - 0).
 *
 * @param s A set of integers of any size with values ranging from 1..UINT_MAX-1
 * @return maximum range between the numbers or 0 and the smallest number
 */
unsigned int MultiEDSM::getMaxRangeInSet(vector<unsigned int> s)
{
    if (s.size() == 1)
    {
        return s[0];
    }
    else if (s.size() == 2)
    {
        unsigned int m = (s[0] > s[1]) ? s[0] - s[1] : s[1] - s[0];
        return max(m, min(s[0], s[1]));
    }
    else if (s.size() == 3)
    {
        unsigned int x[] = {s[0], s[1], s[2]};
        if (s[1] < s[0]) {
            x[0] = s[1];
            x[1] = s[0];
        }
        if (s[2] < x[1]) {
            x[2] = x[1];
            x[1] = s[2];
            if (x[1] < x[0]) {
                unsigned int temp = x[0];
                x[0] = x[1];
                x[1] = temp;
            }
        }
        unsigned int m = max(x[1] - x[0], x[2] - x[1]);
        return max(m, x[0]);
    }
    else
    {
        unsigned int i, n = s.size(), mn = UINT_MAX, mx = 0, MINVAL = 0, MAXVAL = UINT_MAX;
        for (i = 0; i < n; i++) {
            if (s[i] < mn) {
                mn = s[i];
            }
            if (s[i] > mx) {
                mx = s[i];
            }
        }

        unsigned int * minima_bins = new unsigned int[n - 1];
        unsigned int * maxima_bins = new unsigned int[n - 1];
        for (i = 0; i < n - 1; i++) {
            minima_bins[i] = MAXVAL;
            maxima_bins[i] = MINVAL;
        }

        unsigned int bin_idx;
        double bin_size = (double)(mx - 1) / (double)(n - 1);

        for (i = 0; i < n; i++)
        {
            if (s[i] == mn || s[i] == mx) {
                continue;
            }
            bin_idx = (unsigned int) ((double)(s[i] - mn) / bin_size);
            if (minima_bins[bin_idx] == MAXVAL || s[i] < minima_bins[bin_idx]) {
                minima_bins[bin_idx] = s[i];
            }
            if (maxima_bins[bin_idx] == MINVAL || s[i] > maxima_bins[bin_idx]) {
                maxima_bins[bin_idx] = s[i];
            }
        }

        unsigned int prev = mn, max_gap = 0;
        for (i = 0; i < n - 1; i++)
        {
            if (minima_bins[i] == MAXVAL) {
                continue;
            }
            if ((minima_bins[i] - prev) > max_gap) {
                max_gap = minima_bins[i] - prev;
            }
            prev = maxima_bins[i];
        }

        delete [] minima_bins;
        delete [] maxima_bins;

        if ((mx - prev) > max_gap) {
            max_gap = mx - prev;
        }
        if (mn > max_gap) {
            max_gap = mn;
        }

        return max_gap;
    }
}

void MultiEDSM::constructOV7(const string & p)
{
    //
    // Step 1: Traverse the tree level order down to maxP-2 (inclusive) levels
    // making sure not to follow edges beginning with a separator char
    //
    // cout << "Creating nodeLevels..." << endl;
    vector<vector<unsigned int>> nodeLevels;
    vector<unsigned int> v;
    v.push_back(this->STp.id(this->STp.root()));
    nodeLevels.push_back(v);
    queue<cst_node_t> q;
    q.push(this->STp.root());
    unsigned int maxDepth = (unsigned int) max((int)1, (int)this->maxP - 2); //max level to traverse in level order
    unsigned int cumulativeNodeTotal = 0, i = 0; //CNT is maximum nodes to store in either OVMemU7, OVMem7 or between them
    for (i = 1; i <= maxDepth; i++) {
        cumulativeNodeTotal += pow(SIGMA, i) + 1;
    }
    this->OVMemU7.reserve(min((unsigned int)this->STp.nodes(), min(this->maxK, cumulativeNodeTotal))); //maximum number of nodes to store in OVMemU7
    unsigned int numNodesInOVMemU7 = 0; //count of nodes actually stored in OVMemU7
    unsigned int level = 1;             //current level -- root is level 0, but we are inserting children into q, so level = currNodeLevel+1
    unsigned int dc = 1;                //decrement counter
    unsigned int ic = 0;                //increment counter
    unsigned int nn = 0;                //num nodes in nodeLevels
    unsigned int id;                    //node id
    cst_node_t currNode;
    while (!q.empty() && level <= maxDepth)
    {
        currNode = q.front();
        q.pop();
        for (const auto & child : this->STp.children(currNode))
        {
            if (this->STp.is_leaf(child))
            {
                id = this->STp.id(child);
                if (nodeLevels.size() <= level) {
                    vector<unsigned int> v;
                    v.push_back(id);
                    nodeLevels.push_back(v);
                } else {
                    nodeLevels[level].push_back(id);
                }
                if (numNodesInOVMemU7 <= this->maxK) {
                    WordVector w;
                    this->OVMemU7[id] = w;
                    numNodesInOVMemU7++;
                }
                q.push(child);
                nn++;
                ic++;
            }
            else
            {
                unsigned int lb = this->STp.lb(child);
                unsigned int sn = this->STp.csa[lb];
                // char nextChar = p[sn + level - 1];
                char nextChar = p[sn];
                if (nextChar != SEPARATOR_CHAR) {
                    id = this->STp.id(child);
                    if (nodeLevels.size() <= level) {
                        vector<unsigned int> v;
                        v.push_back(id);
                        nodeLevels.push_back(v);
                    } else {
                        nodeLevels[level].push_back(id);
                    }
                    if (numNodesInOVMemU7 <= this->maxK) {
                        WordVector w;
                        this->OVMemU7[id] = w;
                        numNodesInOVMemU7++;
                    }
                    q.push(child);
                    nn++;
                    ic++;
                }
            }
        }
        dc--;
        if (dc == 0) {
            level++;
            dc = ic;
            ic = 0;
        }
    }

    // for (unsigned int x = 0; x < nodeLevels.size(); x++)
    // {
    //     cout << "Level " << x << " contains " << nodeLevels[x].size() << " nodes;" << endl;
    //     for (unsigned int y : nodeLevels[x]) {
    //         cout << y << " ";
    //     }
    //     cout << endl;
    // }
    // With separator char checking (buggy nextChar = p[sn + level - 1]):
    // Level 0 contains 1 nodes.
    // Level 1 contains 5 nodes.
    // Level 2 contains 17 nodes.
    // Level 3 contains 65 nodes.
    // Level 4 contains 257 nodes.
    // Level 5 contains 1025 nodes.
    // Level 6 contains 4097 nodes.
    // Level 7 contains 16385 nodes.
    // Level 8 contains 65537 nodes.
    // Level 9 contains 260952 nodes.
    // With seperator char checking (p[sn]):
    // Level 0 contains 1 nodes.
    // Level 1 contains 5 nodes.
    // Level 2 contains 21 nodes.
    // Level 3 contains 97 nodes.
    // Level 4 contains 449 nodes.
    // Level 5 contains 2049 nodes.
    // Level 6 contains 9217 nodes.
    // Level 7 contains 40961 nodes.
    // Level 8 contains 180225 nodes.
    // Level 9 contains 784202 nodes.
    // Without separator char checking:
    // Level 0 contains 1 nodes.
    // Level 1 contains 6 nodes.
    // Level 2 contains 25 nodes.
    // Level 3 contains 113 nodes.
    // Level 4 contains 513 nodes.
    // Level 5 contains 2305 nodes.
    // Level 6 contains 10241 nodes.
    // Level 7 contains 45057 nodes.
    // Level 8 contains 196609 nodes.
    // Level 9 contains 849631 nodes.


    //
    // Step 2: For each node at maxDepth, encode.
    //
    // cout << "Starting node encoding..." << endl;
    this->OVMem7.reserve(max(1, (int)nn - (int)numNodesInOVMemU7));
    maxDepth = nodeLevels.size() - 1;
    for (unsigned int nodeId : nodeLevels[maxDepth])
    {
        currNode = this->STp.inv_id(nodeId);
        this->OV7EncodeNodes(currNode, nodeId);
    }

    //
    // Step 3: Traverse the tree in reverse level order to encode children into
    // their parents from maxDepth down to level 1.
    //
    // cout << "Tree climbing..." << endl;
    bool currNodeIsInOVMemU7, parentIsInOVMemU7;
    unsigned int parentId, sn, j;
    for (i = maxDepth; i > 0; i--)
    {
        for (unsigned int nodeId : nodeLevels[i])
        {
            currNode = this->STp.inv_id(nodeId);
            currNodeIsInOVMemU7 = (this->OVMemU7.find(nodeId) != this->OVMemU7.end());
            parentId = this->STp.id(this->STp.parent(currNode));
            parentIsInOVMemU7 = (this->OVMemU7.find(parentId) != this->OVMemU7.end());

            if (!currNodeIsInOVMemU7 && !parentIsInOVMemU7)
            {
                // cout << "Both not WordVectors ";
                for (unsigned int x : this->OVMem7[nodeId]) {
                    this->OVMem7[parentId].push_back(x);
                }
                // cout << "nodeVal:" << this->OVMemU7[nodeId] << " parentId:" << parentId << endl;
            }
            else if (!currNodeIsInOVMemU7 && parentIsInOVMemU7)
            {
                // cout << "Only parent is WordVector ";
                for (unsigned int x : this->OVMem7[nodeId]) {
                    this->WordVectorSet1At(this->OVMemU7[parentId], x);
                }
                // cout << "nodeVal:" << this->OVMemU7[nodeId] << " parentId:" << parentId << endl;
            }
            else if (currNodeIsInOVMemU7 && parentIsInOVMemU7)
            {
                if (this->STp.is_leaf(currNode))
                {
                    // cout << nodeId << " Leaf in nodeLevels ";
                    sn = this->STp.sn(currNode);
                    if (
                        (sn > 0 && ((sn + 1) < (this->R - 1))) && \
                        (this->STpIdx2BVIdx[sn-1] != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn]   != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn+1] != SEPARATOR_DIGIT)
                    ) {
                        //calculate bit position
                        j = this->STpIdx2BVIdx[sn] - 1;  //j is index in bitvector
                        this->WordVectorSet1At(this->OVMemU7[nodeId], j);
                        // cout << "nodeVal:" << this->OVMemU7[nodeId] << " parentId:" << parentId << endl;
                    }
                    else
                    {
                        // cout << "not encoded, sn:" << sn << endl;
                    }
                }
                else
                {
                    // cout << nodeId << " Not leaf ";
                }
                this->WordVectorOR_IP(this->OVMemU7[parentId], this->OVMemU7[nodeId]);
                // cout << "nodeVal:" << this->OVMemU7[nodeId] << " parentId:" << parentId << endl;
            }
            else if (!parentIsInOVMemU7 && this->STp.is_leaf(currNode))
            {
                if (this->STp.is_leaf(currNode))
                {
                    // cout << nodeId << " Leaf in nodeLevels but not parent ";
                    sn = this->STp.sn(currNode);
                    if (
                        (sn > 0 && ((sn + 1) < (this->R - 1))) && \
                        (this->STpIdx2BVIdx[sn-1] != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn]   != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn+1] != SEPARATOR_DIGIT)
                    ) {
                        //calculate bit position
                        j = this->STpIdx2BVIdx[sn] - 1;  //j is index in bitvector
                        this->WordVectorSet1At(this->OVMemU7[nodeId], j);
                        // cout << "nodeVal:" << this->OVMemU7[nodeId] << " parentId:" << parentId << endl;
                    }
                    else
                    {
                        // cout << "not encoded, nodeId:" << nodeId << endl;
                    }
                }
                else
                {
                    // cout << nodeId << " Not leaf ";
                }
            }
            else
            {
                //Node not useful so ignored. Or currently at level 1, where
                //parent is root and children are WordVectors; nothing to do here.
                // cout << nodeId << " Nope ";
                // cout << "nodeVal:" << this->OVMemU7[nodeId] << " parentId:" << parentId << endl;
            }
        }
        nodeLevels[i].clear();
    }
    nodeLevels.clear();

    // cout << endl << "idx\tsn\tsuffix" << endl;
    // for (i = 0; i < p.length() + 1; i++)
    // {
    //     cout << i << "\t" << this->STp.csa[i] << "\t" << p.substr(this->STp.csa[i]) << endl;
    // }
    // cout << endl;
}
void MultiEDSM::OV7EncodeNodes(const cst_node_t & currNode, const unsigned int nodeId)
{
    if (this->STp.is_leaf(currNode))
    {
        // cout << "Encoding leaf " << nodeId << " ";
        unsigned int sn = this->STp.sn(currNode);
        unsigned int j = this->STpIdx2BVIdx[sn] - 1;
        if (this->OVMemU7.find(nodeId) != this->OVMemU7.end() && sn > 0) {
            this->WordVectorSet1At(this->OVMemU7[nodeId], j);
            this->OVMemU7[nodeId].shrink_to_fit();
            // cout << "sn" << sn << " pos" << j << " " << this->OVMemU7[nodeId];
        } else {
            if (sn > 0) {
                vector<unsigned int> v;
                v.push_back(j);
                v.shrink_to_fit();
                this->OVMem7[nodeId] = v;
                // cout << v;
            }
        }
        // cout << endl;
    }
    else
    {
        // cout << "Encoding node " << nodeId << " ";
        unsigned int lb = this->STp.lb(currNode);
        unsigned int rb = this->STp.rb(currNode);
        unsigned int sn, i, j;
        bool isInOVMemU7 = (this->OVMemU7.find(nodeId) != this->OVMemU7.end());
        for (i = lb; i < rb + 1; i++)
        {
            sn = this->STp.csa[i];
            if (
                (sn > 0 && ((sn + 1) < (this->R - 1))) && \
                (this->STpIdx2BVIdx[sn-1] != SEPARATOR_DIGIT) && \
                (this->STpIdx2BVIdx[sn]   != SEPARATOR_DIGIT) && \
                (this->STpIdx2BVIdx[sn+1] != SEPARATOR_DIGIT)
            ) {
                //calculate bit position
                j = this->STpIdx2BVIdx[sn] - 1;  //j is index in bitvector
                if (isInOVMemU7) {
                    this->WordVectorSet1At(this->OVMemU7[nodeId], j);
                    // cout << this->OVMemU7[nodeId] << endl;
                } else {
                    this->OVMem7[nodeId].push_back(j);
                    // cout << this->OVMem7[nodeId] << endl;
                }
            }
        }
        if (!isInOVMemU7) {
            this->OVMemU7[nodeId].shrink_to_fit();
        } else {
            // cout << "Size of Node " << nodeId << " is " << (this->OVMem7[nodeId].size() * sizeof(unsigned int)) << endl;
            this->OVMem7[nodeId].shrink_to_fit();
        }
    }
}


/**
 * The Shift-And algorithm performs the equivalent of creating a BorderPrefixTable
 * and this function just returns a WordVector with the prefixes of patterns
 * identified in the suffixes of all strings in a segment.
 * If the total length of all patterns is M, this is encoded into O([M/w]) computer
 * words. The total length of the strings in a segment is N, so the building of
 * the prefix border table should take time O(N[M/w]). But we need to factor in
 * the fact we do a bitwise OR operation for j strings in a segment, so the time
 * taken is worst case O(N*[jM/w]).
 *
 * @param S A segment
 * @return A WordVector with the prefixes marked
 */
WordVector MultiEDSM::buildBorderPrefixWordVector(const Segment & S)
{
    WordVector c;
    unsigned int m, i = 0;
    Segment::const_iterator stringI;

    this->umsa->disableReporting();

    for (stringI = S.begin(); stringI != S.end(); ++stringI)
    {
        if (stringI[0] != EPSILON)
        {
            m = (*stringI).length();
            if (m < this->maxP) {
                this->umsa->search(*stringI);
            } else {
                this->umsa->search(*stringI, m - this->maxP + 1);
            }
            if (i == 0) {
                c = this->umsa->getLastSearchState();
            } else {
                this->WordVectorOR_IP(c, this->umsa->getLastSearchState());
            }
            i++;
        }
    }

    this->umsa->enableReporting();

    return c;
}

/**
 * occVector tool returns the index in MultiEDSM::OVMem of the starting positions
 * of a given substring _a_, if it is present in the pattern, encoded as a bitvector.
 * If _a_ is not present in any patterns as an infix then 0 will be returned (0
 * index is the root element but it is not used directly by the algorithm).
 * Time taken: O(|a|).
 *
 * @param a A substring
 * @param B2 A WordVector to AND with
 * @return Starting positions of substring _a_ encoded into a bitvector
 */
void MultiEDSM::occVector(const string & a, WordVector & B2)
{
    cst_node_t explicitNode = this->STp.root();
    string::const_iterator it;
    uint64_t char_pos = 0;
    unsigned int j = 0;
    for (it = a.begin(); it != a.end(); ++it)
    {
        if (forward_search(this->STp, explicitNode, it - a.begin(), *it, char_pos) > 0) {
            j++;
        } else {
            break;
        }
    }

    if (j == a.length())
    {
        // cout << a << endl;
        unsigned int nodeId = this->STp.id(explicitNode);
        if (this->OVMemU7.find(nodeId) != this->OVMemU7.end())
        {
            // cout << a << "; B2 before:" << B2 << " ";
            this->WordVectorAND_IP(B2, this->OVMemU7[nodeId]);
            // cout << "OVMemU7[" << nodeId << "]:" << this->OVMemU7[nodeId] << " " \
            //      << "B2 After:" << B2 << endl;
        }
        else if (this->OVMem7.find(nodeId) != this->OVMem7.end())
        {
            WordVector v;
            for (unsigned int pos : this->OVMem7[nodeId]) {
                this->WordVectorSet1At(v, pos);
            }
            this->WordVectorAND_IP(B2, v);
        }
        else
        {
            B2.clear();
        }
    }
    else
    {
        B2.clear();
    }
}

/**
 * Report match found. The results can be obtained by calling MultiEDSM::getMatches().
 * Note: Although three arguments are taken by this method, we only report two
 * of them: the Position and the PatternId. The posIdx parameter, the relative
 * position inside a segment may be useful for some purposes but we do not
 * currently report it.
 *
 * @param matchIdx The index of the ending position in the EDT where the match was found
 * @param posIdx The position inside a segment where the match was found
 * @param pattId The id of the pattern found
 */
void MultiEDSM::report(const unsigned int matchIdx, const unsigned int posIdx, const int pattId)
{
    if (this->matches.size() == 0)
    {
        this->matches.push_back(pair<unsigned int, unsigned int>(matchIdx, pattId));
    }
    else
    {
        pair<unsigned int, unsigned int> & last = this->matches.back();
        if (last.first != matchIdx || last.second != pattId) {
            this->matches.push_back(pair<unsigned int, unsigned int>(matchIdx, pattId));
        }
    }
}

/**
 * Search for the patterns in a segment
 * Time and space complexity:
 *  - O(N*[M/w]) time for long Shift-And search
 *  - O(N*[M/w]) time for BorderPrefixTable construction
 *  - O([M/w]) time for epsilon
 *  - O(N*[M/w]) worst-case time for suffix matching in step 2
 *  - O([M/w]) worst-case time for step 3 bitwise-OR operation and O(min((M + [M/w]), (M([M/w]))) for LeftShift
 *  - O([M/w]) extra space
 */
bool MultiEDSM::searchNextSegment(const Segment & S)
{
    //start timer
    clock_t start = clock();

    //match variables
    bool matchFound = false;
    bool isDeterminateSegment = false;
    unsigned int matchIdx, posIdx;
    int pattId;

    //temp state variables, segment string iterator
    unsigned int k, m;
    WordVector B1, B2;
    Segment::const_iterator stringI;

    // cout << "Position #" << (this->d + this->D) << endl;

    //search first segment, priming B, then search the next segments normally in the else condition down below
    if (!this->primed)
    {
        cout << "Starting search..." << endl;

        //keep track of f/F counts and if segment is determinate
        if (S.size() == 1)
        {
            this->f += S[0].length();
            isDeterminateSegment = true;
        }
        else
        {
            for (stringI = S.begin(); stringI != S.end(); ++stringI)
            {
                if (stringI[0] != EPSILON) {
                    this->F += (*stringI).length();
                }
            }
            isDeterminateSegment = false;
        }

        //search for the patterns
        for (stringI = S.begin(); stringI != S.end(); ++stringI)
        {
            if (stringI[0] != EPSILON)
            {
                if ((*stringI).length() >= this->minP)
                {
                    if (this->umsa->search(*stringI))
                    {
                        if (isDeterminateSegment)
                        {
                            for (const pair<int,int> & match : this->umsa->getMatches()) {
                                posIdx = match.first;
                                matchIdx = this->pos + posIdx;
                                pattId = match.second;
                                this->report(matchIdx, posIdx, pattId);
                            }
                            this->umsa->clearMatches();
                        }
                        else
                        {
                            if (this->reportOnce && !matchFound)
                            {
                                const pair<int,int> match = *this->umsa->getMatches().begin();
                                posIdx = match.first;
                                matchIdx = this->pos;
                                pattId = match.second;
                                this->report(matchIdx, posIdx, pattId);
                                this->umsa->clearMatches();
                                break;
                            }
                            else if (!this->reportOnce)
                            {
                                for (const pair<int,int> & match : this->umsa->getMatches()) {
                                    posIdx = match.first;
                                    matchIdx = this->pos;
                                    pattId = match.second;
                                    this->report(matchIdx, posIdx, pattId);
                                }
                                this->umsa->clearMatches();
                            }
                        }
                        matchFound = true;
                    }
                }
            }
        }

        //build equivalent of BorderPrefixTable for next segment to use
        this->B = this->buildBorderPrefixWordVector(S);
        // cout << "B:" << this->B << endl;

        //update position
        if (isDeterminateSegment) {
            this->d++;
            this->pos = this->f;
        } else {
            this->D++;
            this->pos++;
        }

        this->primed = true;
    }
    else
    {
        //build equivalent of BorderPrefixTable for current segment
        // cout << "B: " << this->B << endl;
        B1 = this->buildBorderPrefixWordVector(S);
        // cout << "B1: " << B1 << endl;

        //keep track of f/F counts and if segment is determinate, and if there's an epsilon (deletion)
        if (S.size() == 1)
        {
            this->f += S[0].length();
            isDeterminateSegment = true;
        }
        else
        {
            for (stringI = S.begin(); stringI != S.end(); ++stringI)
            {
                if (stringI[0] == EPSILON) {
                    this->WordVectorOR_IP(B1, this->B);
                } else {
                    this->F += (*stringI).length();
                }
            }
            isDeterminateSegment = false;
        }

        //loop through the strings in the segment
        for (stringI = S.begin(); stringI != S.end(); ++stringI)
        {
            if (stringI[0] != EPSILON)
            {
                m = (*stringI).length();

                //step 1 - if string is longer than pattern then do a full search on it
                if (m >= this->minP)
                {
                    if (this->umsa->search(*stringI))
                    {
                        if (isDeterminateSegment)
                        {
                            for (const pair<int,int> & match : this->umsa->getMatches()) {
                                posIdx = match.first;
                                matchIdx = this->pos + posIdx;
                                pattId = match.second;
                                this->report(matchIdx, posIdx, pattId);
                            }
                            this->umsa->clearMatches();
                        }
                        else
                        {
                            if (this->reportOnce && !matchFound)
                            {
                                const pair<int,int> match = *this->umsa->getMatches().begin();
                                posIdx = match.first;
                                matchIdx = this->pos;
                                pattId = match.second;
                                this->report(matchIdx, posIdx, pattId);
                                this->umsa->clearMatches();
                            }
                            else if (!this->reportOnce)
                            {
                                for (const pair<int,int> & match : this->umsa->getMatches()) {
                                    posIdx = match.first;
                                    matchIdx = this->pos;
                                    pattId = match.second;
                                    this->report(matchIdx, posIdx, pattId);
                                }
                                this->umsa->clearMatches();
                            }
                        }
                        matchFound = true;
                    }
                }

                //step 2 - if string is a suffix of a previously determined prefix or infix, then report a match
                B2 = this->B;
                if (this->umsa->search(*stringI, B2, 0, this->maxP - 1))
                {
                    if (this->reportOnce && !matchFound)
                    {
                        const pair<int,int> match = *this->umsa->getMatches().begin();
                        posIdx = match.first;
                        if (isDeterminateSegment) {
                            matchIdx = this->pos + posIdx;
                        } else {
                            matchIdx = this->pos;
                        }
                        pattId = match.second;
                        this->report(matchIdx, posIdx, pattId);
                        this->umsa->clearMatches();
                    }
                    else if (!this->reportOnce)
                    {
                        for (const pair<int,int> & match : this->umsa->getMatches()) {
                            posIdx = match.first;
                            if (isDeterminateSegment) {
                                matchIdx = this->pos + posIdx;
                            } else {
                                matchIdx = this->pos;
                            }
                            pattId = match.second;
                            this->report(matchIdx, posIdx, pattId);
                        }
                        this->umsa->clearMatches();
                    }
                    matchFound = true;
                }

                //step 3 - if string is short, then it could be an infix, so determine it is actual infix and mark its position accordingly
                if (m < (this->maxP - 1))
                {
                    this->Np++;
                    this->Nm += m;
                    B2 = this->B;
                    this->occVector(*stringI, B2);
                    if (B2.size() > 0) {
                        this->WordVectorLeftShift_IP(B2, m);
                        this->WordVectorOR_IP(B1, B2);
                    }
                }
            }
        }

        //update the state for the next segment search
        this->B = B1;

        //update position
        if (isDeterminateSegment) {
            this->d++;
            this->pos += S[0].length();
        } else {
            this->D++;
            this->pos++;
        }
    }

    //stop timer
    this->duration += clock() - start;

    return matchFound;
}

/**
 * The bitwise OR operation performed 'in place' on parameter a.
 *
 * @param a Subject WordVector which will be altered
 * @param b Object WordVector
 */
void MultiEDSM::WordVectorOR_IP(WordVector & a, const WordVector & b)
{
    unsigned int i, as = a.size(), bs = b.size();
    if (as == 0)
    {
        a = b;
    }
    else if (as == bs) // && as != 0
    {
        for (i = 0; i < as; i++) {
            a[i] = a[i] | b[i];
        }
    }
    else if (as > bs)
    {
        for (i = 0; i < bs; i++) {
            a[i] = a[i] | b[i];
        }
    }
    else
    {
        for (i = 0; i < as; i++) {
            a[i] = a[i] | b[i];
        }
        for (; i < bs; i++) {
            a.push_back(b[i]);
        }
    }
}

/**
 * The bitwise AND operation performed 'in place' on parameter a. Time taken:
 * O([M/w])
 *
 * @param a WordVector subject which will be altered
 * @param b WordVector object
 */
void MultiEDSM::WordVectorAND_IP(WordVector & a, const WordVector & b)
{
    unsigned int i, m = min(a.size(), b.size()), n = a.size();
    for (i = 0; i < m; i++) {
        a[i] = a[i] & b[i];
    }
    for (; i < n; i++) {
        a[i] = 0;
    }
}

/**
 * A special left shifting algorithm for use with OccVector result. It ensures a
 * shift operation is not performed if it overlaps into the next pattern in the
 * bitvector. This method uses the Pos2PatId datastructure to identify if a shift
 * illegally crosses into the next pattern, so takes O(M) space.
 * Time taken is roughly: O(num_set_bits(x) + [M/w])... or really
 * worst-Worst-WORST case, where there is a match at every single position, which
 * is very-Very-VERY unlikely then O(M + [M/w])!
 *
 * @param x
 * @param m The length of the string being checked
 * @TODO figure out why this method is buggy and restore in MultiEDSM::WordVectorLeftShift_IP
 */
void MultiEDSM::WordVectorSPECIALSHIFT_IP(WordVector & x, unsigned int m)
{
    WORD temp;
    unsigned int i, j, currPos, updPos, currPattId, updPattId, currWordIdx, updWordIdx, n = x.size();
    for (i = 0; i < n; i++)
    {
        temp = x[i];
        while (temp)
        {
            j = ffs(temp) - 1;
            currPos = i * BITSINWORD + j;
            currPattId = this->Pos2PatId[currPos];
            updPos = currPos + m;
            updPattId = this->Pos2PatId[updPos];
            if (currPattId == updPattId)
            {
                //remove old bit position
                // currWordIdx = (unsigned int) ((float)currPos / (float)BITSINWORD);
                // x[currWordIdx] = x[currWordIdx] ^ (1ul << j);
                //set new bit position
                updWordIdx = (unsigned int) ((double)updPos / (double)BITSINWORD);
                while (updWordIdx >= n) {
                    x.push_back(0ul);
                }
                x[updWordIdx] = x[updWordIdx] | (1ul << (updPos % BITSINWORD));
            }
            temp = temp ^ (1ul << j);
        }
    }
}

/**
 * A simple left shifting algorithm for use with OccVector result. It shifts all
 * the words in a WordVector m times, making sure the moving 1's don't overspill
 * past the end position of a pattern into the next pattern. Time taken is
 * O(m[M/w]), and does not use up any extra space because the ending positions
 * are obtained from UMSA's Shift-And algorithm, occupying O([M/w]) space.
 *
 * @param x
 * @param m The length of the string being checked
 */
void MultiEDSM::WordVectorSIMPLESHIFT_IP(WordVector & x, unsigned int m)
{
    unsigned int k = x.size();

    if (k > 0)
    {
        const WordVector & ends = this->umsa->getEndingStates();
        unsigned int i, j;
        WORD temp, carry, carryMask = 1ul << (BITSINWORD - 1);
        for (i = 0; i < m; i++)
        {
            carry = 0;
            for (j = 0; j < k; j++)
            {
                temp = x[j];
                x[j] = (x[j] << 1) | carry; //left shift and apply the carried 1
                x[j] = x[j] ^ (ends[j] & x[j]); //if 1 passes pattern boundary (end), then remove it because it is an illegal shift
                carry = (WORD)((carryMask & temp) != 0); //work out if we need to carry a 1 to the next word
                if (carry && (j == (k - 1))) { //identify if we need to expand x //@TODO maybe this is not required because x will always be the correct size?
                   x.push_back(0ul);
                   k++;
                }
            }
        }
    }
}

/**
 * Performs the left-shift operation for the WordVector returned after running the
 * OccVector tool. This operation is potentially expensive so there are two
 * different ways to do it:
 *  - MultiEDSM::WordVectorSPECIALSHIFT_IP() has a worst case time of O(num_set_bits(x) + [M/w])
 *    or O(M + [M/w]) if there is a huge number of set bits in x.
 *  - MultiEDSM::WordVectorSIMPLESHIFT_IP() has a worst case time of O(m[M/w])
 *    because it does m shift operations over [M/w] words.
 * To get the best performance, we do a triage and identify which of the methods
 * does it in the least time. The triage check itself takes O([M/w]) time.
 *
 * @param x
 * @param m The length of the string being checked
 */
void MultiEDSM::WordVectorLeftShift_IP(WordVector & x, unsigned int m)
{
    // unsigned int i = 0, p = 0, s = x.size();
    // unsigned int specialShiftTime = 0, simpleShiftTime = m * s;
    //
    // while (i < s && specialShiftTime <= simpleShiftTime)
    // {
    //     p += popcount(x[i]);
    //     specialShiftTime = p + s;
    //     i++;
    // }

    // if (simpleShiftTime <= specialShiftTime) {
        this->WordVectorSIMPLESHIFT_IP(x, m);
    // } else {
    //     this->WordVectorSPECIALSHIFT_IP(x, m);
    // }
}

/**
 * Sets the bit at position pos to 1 in WordVector x. 0-based bit index. If x is
 * too small it will be expanded to accomodate the new bit.
 *
 * @param x WordVector
 * @param pos The position to place the 1 at
 */
void MultiEDSM::WordVectorSet1At(WordVector & x, unsigned int pos)
{
    unsigned int wordIdx = (unsigned int) ((double)pos / (double)BITSINWORD);
    unsigned int bitShifts = pos % BITSINWORD;
    unsigned int j = x.size();

    while (j <= wordIdx) {
        x.push_back(0ul);
        j++;
    }

    x[wordIdx] = x[wordIdx] | (1ul << bitShifts);
}

/**
* Get a list of the positions of matches found.
*
* @return Matches found
*/
ResultSet MultiEDSM::getMatches() const
{
    return this->matches;
}

/**
* Clear the matches list
*/
void MultiEDSM::clearMatches()
{
    this->matches.clear();
}

/**
 * Report the position once only in indeterminate segments?
 *
 * @param yesorno Default: True
 */
void MultiEDSM::reportOncePerPosition(bool yesorno)
{
    this->reportOnce = yesorno;
}

/**
 * Report the id of discovered patterns?
 *
 * @param yesorno Default: True
 */
void MultiEDSM::reportPatternIds(bool yesorno)
{
    this->reportPatterns = yesorno;
}

/**
* Get the total number of determinate segments searched so far
*/
unsigned int MultiEDSM::getd() const
{
    return this->d;
}

/**
* Get the total number of degenerate segments searched so far
*/
unsigned int MultiEDSM::getD() const
{
    return this->D;
}

/**
* Get the total length of determinate segments searched so far
*/
unsigned int MultiEDSM::getf() const
{
    return this->f;
}

/**
* Get the total length of degenerate segments searched so far
*/
unsigned int MultiEDSM::getF() const
{
    return this->F;
}

/**
* Get the total number of strings analyzed that are shorter than minP
*/
unsigned int MultiEDSM::getNp() const
{
    return this->Np;
}

/**
* Get the total length of the strings analyzed that are shorter than minP
*/
unsigned int MultiEDSM::getNm() const
{
    return this->Nm;
}

/**
* Returns the execution duration of MultiEDSM in seconds
*/
double MultiEDSM::getDuration() const
{
    return this->duration / (double) CLOCKS_PER_SEC;
}
