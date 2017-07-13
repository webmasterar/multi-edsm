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
 */
MultiEDSM::MultiEDSM(const string & alphabet, const vector<string> & patterns)
{
    this->alphabet = alphabet;
    if (patterns.size() == 0) {
        cerr << "Error: Empty pattern set!" << endl;
        throw 1;
    }
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
    char sep = '#';
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

        //add pattern to p and append with seperator for building suffix tree of patterns
        p += patterns[i];
        p += sep;

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

    //construct the suffix tree of patterns with seperators
    cout << "SuffixTree..." << endl;
    construct_im(this->STp, p.c_str(), sizeof(char));

    //tool to get pattern id from any given position in the bitvector/p
    cout << "Pos2Pat" << endl;
    this->Pos2PatId = this->umsa->getPatternPositions();

    //construct occVector datastructure
    this->constructOV4();

    //stop timer
    this->duration += clock() - start;

    // cout << "Preproccessing " << patterns.size() << " patterns of size M=" \
    //      << this->M << " completed in " << this->getDuration() << "s." << endl;
}

void MultiEDSM::constructOV4()
{
    unsigned int i, numNodes = this->STp.nodes();
    cout << "NumNodes: " << numNodes << endl;
    for (i = 0; i < numNodes; i++) {
        WordVector v(1, 0ul);
        v.shrink_to_fit();
        this->OVMem.push_back(v);
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
            this->OVMem[id] = a;
        }
    }
    else
    {
        for (const auto & child : this->STp.children(u)) {
            this->WordVectorOR_IP(this->OVMem[id], this->recAssignOVMem4(child));
        }
    }

    return this->OVMem[id];
}

/**
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

                    //make sure it is not first letter or last letter in a pattern
                    if (h > 0 && this->STpIdx2BVIdx[h - 1] != SEPARATOR_DIGIT && \
                       (h + 1) < ((int)this->R - 1) && this->STpIdx2BVIdx[h + 1] != SEPARATOR_DIGIT)
                    {
                        //calculate bit position
                        long long int k = this->STpIdx2BVIdx[h] - 1;  //i is index in bitvector
                        unsigned int l = 1;                           //current size of a
                        unsigned int wordIdx = (unsigned int) ((double)k / (double)BITSINWORD);
                        unsigned int bitShifts = k % BITSINWORD;

                        //create word vector with enough words in it
                        WordVector a(wordIdx + 1, 0ul);

                        //write position to a, OVMem[id] if applicable and OVMem[parent_id]
                        a[wordIdx] = a[wordIdx] | (1ul << bitShifts);
                        if (j <= this->maxP) {
                            this->OVMem2[id] = a;
                        }

                        //create the parent if not exists or combine a with parent
                        if (this->OVMem2.find(parent_id) == this->OVMem2.end()) {
                            this->OVMem2[parent_id] = a;
                            leafCreationCounter++;
                            totalLeafSize += wordIdx + 1;
                        } else {
                            this->WordVectorOR_IP(this->OVMem2[parent_id], a);
                        }
                    }
                }
                else
                {
                    //if node already exists
                    if (this->OVMem2.find(id) != this->OVMem2.end())
                    {
                        if (this->OVMem2.find(parent_id) == this->OVMem2.end()) {
                            this->OVMem2[parent_id] = this->OVMem2[id];
                            nodeCreationCounter++;
                        } else {
                            this->WordVectorOR_IP(this->OVMem2[parent_id], this->OVMem2[id]);
                        }
                        if (this->STp.depth(currNode) > this->maxP) {
                            // this->OVMem[id].empty();
                            this->OVMem2.erase(id);
                        }
                    }
                }
            }
        }

        totalNodesCounter = nodeLevels[j].size();

        cout << "Level " << j << " contains " << totalNodesCounter << " nodes. " \
             << leafCreationCounter << " leaves and " << nodeCreationCounter \
             << " nodes were created. Average leaf size in bytes: " \
             << (8 * (int)((double)totalLeafSize / (double)leafCreationCounter)) << endl;

        nodeLevels[j].empty();
    }
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
 * @TODO rewrite
 */
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
 * @return Starting positions of substring _a_ encoded into a bitvector
 */
unsigned int MultiEDSM::occVector(const string & a)
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

    return (j == a.length()) ? this->STp.id(explicitNode) : 0;
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

    //search first segment, priming B, then search the next segments normally in the else condition
    if (!this->primed)
    {
        cout << "Starting search" << endl;

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
        B1 = this->buildBorderPrefixWordVector(S);

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
                    k = this->occVector(*stringI);
                    if (k != 0)
                    {
                        B2 = this->B;
                        this->WordVectorAND_IP(B2, this->OVMem[k]);
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
    const WordVector & ends = this->umsa->getEndingStates();

    WORD temp, carry, carryMask = 1ul << (BITSINWORD - 1);
    unsigned int i, j, k = x.size();
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
