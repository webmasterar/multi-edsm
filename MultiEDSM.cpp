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
 *  - OccVector construction time O(M + k), space O((maxK-2) * [(M+k)/w])
 *
 * @param patterns
 */
void MultiEDSM::preprocessPatterns(const vector<string> & patterns)
{
    cout << "Preprocessing starting..." << endl;

    //start timer
    clock_t start = clock();

    //initialize Shift-And searching tool
    cout << "1. Shift-And" << endl;
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

        //add pattern to p and append with separator (#) for building suffix tree of patterns
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
    p.pop_back(); //remove unnecessary terminal separator char as \0 already tagged on end of c_str

    //tool to get pattern id from any given position in the bitvector/p
    cout << "2. Pos2Pat" << endl;
    this->Pos2PatId = this->umsa->getPatternPositions();

    //construct the suffix tree of patterns with separators
    cout << "3. SuffixTree" << endl;
    construct_im(this->STp, p.c_str(), sizeof(char));

    //construct occVector datastructure
    cout << "4. OccVector" << endl;
    this->constructOV8(p);

    //stop timer
    this->duration += clock() - start;

    cout << "Preproccessing " << patterns.size() << " patterns of total size M=" \
         << this->M << " completed in " << this->getDuration() << "s." << endl << endl;
}

/*
* Construct the occVector datastructure. Requires O(M + k) time (i.e. the number
* of nodes in a suffix tree is at most 2M, and k is the number of separators).
* Requires O([(M + k) / w] * maxK) space, where maxK is the maximum number of
* bitvectors to use in the datastructure as determined by what is left of the
* memory specified by the user. It is possible less than maxK explicit nodes
* exist in the suffix tree in which case fewer than maxK bitvectors are used.
*
* @param p The pattern string which the suffix tree was built on of size R = M + k
*/
void MultiEDSM::constructOV8(const string & p)
{
    //
    // Step 1: Traverse the tree level order down to maxP-2 (inclusive) levels
    // keeping a record of which node is found on what level.
    //
    vector<vector<unsigned int>> nodeLevels;
    vector<unsigned int> v;
    v.push_back(this->STp.id(this->STp.root()));
    queue<cst_node_t> q;
    q.push(this->STp.root());
    unsigned int maxDepth = (unsigned int) max((int)1, (int)this->maxP - 2); //max level to traverse in level order
    //Since OVMemU8 is used to store only internal nodes, we subtract R (no. of leaves) from total nodes, to get the max no. If there is insufficient memory only maxK items will be stored
    this->OVMemU8.reserve(min((unsigned int)this->STp.nodes() - this->R, this->maxK)); //maximum number of nodes to store in OVMemU8
    unsigned int maxNoWordsInVector = (unsigned int) ceil((double)this->M / (double)BITSINWORD);
    unsigned int numNodesInOVMemU8 = 0; //count of nodes actually stored in OVMemU8
    unsigned int level = 0;             //current level -- root is level 0
    unsigned int dc = 1;                //decrement counter
    unsigned int ic = 0;                //increment counter
    unsigned int nn = 0;                //num nodes in nodeLevels
    unsigned int id;                    //node id
    cst_node_t currNode;
    bool currNodeIsLeaf;
    while (!q.empty() && level <= maxDepth)
    {
        currNode = q.front();
        currNodeIsLeaf = this->STp.is_leaf(currNode);
        q.pop();

        id = this->STp.id(currNode);
        if (nodeLevels.size() <= level) {
            vector<unsigned int> v;
            v.push_back(id);
            nodeLevels.push_back(v);
        } else {
            nodeLevels[level].push_back(id);
        }
        if (nn > 0 && numNodesInOVMemU8 < this->maxK && !currNodeIsLeaf) { //no leaves allowed in OVMemU8
            WordVector w(maxNoWordsInVector, 0ul);
            this->OVMemU8[id] = w;
            numNodesInOVMemU8++;
        }

        //handle root node special condition -- don't add $ and # starting nodes in first row
        if (nn == 0)
        {
            for (const auto & child : this->STp.children(currNode))
            {
                unsigned int sn, lb;
                if (this->STp.is_leaf(child))
                {
                    sn = this->STp.sn(child);
                }
                else
                {
                    lb = this->STp.lb(child);
                    sn = this->STp.csa[lb];
                }

                if (!(sn == p.length() || p[sn] == SEPARATOR_CHAR))
                {
                    q.push(child);
                    nn++;
                    ic++;
                }
            }
        }
        //handle rest of nodes from level 2 onwards
        else if (!currNodeIsLeaf)
        {
            for (const auto & child : this->STp.children(currNode))
            {
                q.push(child);
                nn++;
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
    //specifying how many nodes/leaves we need to store in OVMem8
    this->OVMem8.reserve(this->STp.nodes() - numNodesInOVMemU8);

    //
    // Step 2: For each explicit node at level maxDepth, encode it, so we can be
    // ready for step 3.
    //
    maxDepth = nodeLevels.size() - 1;
    for (unsigned int nodeId : nodeLevels[maxDepth])
    {
        currNode = this->STp.inv_id(nodeId);
        if (!this->STp.is_leaf(currNode))
        {
            unsigned int lb = this->STp.lb(currNode);
            unsigned int rb = this->STp.rb(currNode);
            if (this->OVMemU8.find(nodeId) != this->OVMemU8.end())
            {
                unsigned int i, j, sn;
                for (i = lb; i <= rb; i++)
                {
                    sn = this->STp.csa[i];
                    if (
                        (sn > 0 && ((sn + 1) < (this->R - 1))) && \
                        (this->STpIdx2BVIdx[sn-1] != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn]   != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn+1] != SEPARATOR_DIGIT)
                    ) {
                        j = this->STpIdx2BVIdx[sn] - 1;  ///calculate bit position -- j is index in bitvector
                        this->WordVectorSet1At(this->OVMemU8[nodeId], j);
                    }
                }
            }
            else
            {
                struct NodeRange nr;
                nr.start = lb;
                nr.end = rb;
                this->OVMem8[nodeId] = nr;
            }
        }
    }

    //
    // Step 3: Climbing up the tree in reverse level order to encode children into
    // their parents from maxDepth down to level 1. This step bitwize ORs explicit
    // nodes with their parent or just specifies their NodeRange and possibly
    // writes the positions to its parent. Or, if the current node is a leaf, it
    // specifies its position and possibly writes the position to its parent.
    //
    bool currNodeIsInOVMemU8, parentIsInOVMemU8;
    unsigned int parentId, sn, i, j, k;
    for (i = maxDepth; i > 0; i--)
    {
        for (unsigned int nodeId : nodeLevels[i])
        {
            currNode = this->STp.inv_id(nodeId);
            currNodeIsInOVMemU8 = (this->OVMemU8.find(nodeId) != this->OVMemU8.end());
            currNodeIsLeaf = this->STp.is_leaf(currNode);
            parentId = this->STp.id(this->STp.parent(currNode));
            parentIsInOVMemU8 = (this->OVMemU8.find(parentId) != this->OVMemU8.end());

            //leaves are assigned in OVMem8
            if (currNodeIsLeaf)
            {
                struct NodeRange nr;
                nr.start = nodeId;
                nr.end = nodeId;
                this->OVMem8[nodeId] = nr;
            }

            //if neither currNode or its parent is in OVMemU8, assign currNode in OVMem8
            if (!currNodeIsInOVMemU8 && !parentIsInOVMemU8 && !currNodeIsLeaf)
            {
                struct NodeRange nr;
                nr.start = this->STp.lb(currNode);
                nr.end = this->STp.rb(currNode);
                this->OVMem8[nodeId] = nr;
            }
            //when both nodes in OVMemU8, just bitwize OR child into parent
            else if (currNodeIsInOVMemU8 && parentIsInOVMemU8)
            {
                this->WordVectorOR_IP(this->OVMemU8[parentId], this->OVMemU8[nodeId]);
            }
            //when only the parent is in OVMemU8, encode child leaves or non-OVMemU8 nodes into parent
            else if (parentIsInOVMemU8 && !currNodeIsInOVMemU8)
            {
                for (j = this->OVMem8[nodeId].start; j <= this->OVMem8[nodeId].end; j++)
                {
                    sn = this->STp.csa[j];
                    if (
                        (sn > 0 && ((sn + 1) < (this->R - 1))) && \
                        (this->STpIdx2BVIdx[sn-1] != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn]   != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn+1] != SEPARATOR_DIGIT)
                    ) {
                        k = this->STpIdx2BVIdx[sn] - 1;  //k is index in bitvector
                        this->WordVectorSet1At(this->OVMemU8[parentId], k);
                    }
                }
            }
        }
        nodeLevels[i].clear();
    }
    nodeLevels.clear();
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
 * @param c A reference to a WordVector to mark prefixes on
 */
void MultiEDSM::buildBorderPrefixWordVector(const Segment & S, WordVector & c)
{
    this->umsa->disableReporting();

    unsigned int m, h, i = 0;
    Segment::const_iterator stringI;
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
                if (c.size() == 0) {
                    c = this->umsa->getLastSearchState();
                } else {
                    const WordVector & s = this->umsa->getLastSearchState();
                    for (h = 0; h < s.size(); h++) {
                        c[h] = s[h];
                    }
                }
            } else {
                this->WordVectorOR_IP(c, this->umsa->getLastSearchState());
            }
            i++;
        }
    }

    this->umsa->enableReporting();
}

/**
 * occVector tool returns the index in MultiEDSM::OVMem of the starting positions
 * of a given substring _a_, if it is present in the pattern, encoded as a bitvector.
 * If _a_ is not present in any patterns as an infix then 0 will be returned (0
 * index is the root element but it is not used directly by the algorithm).
 * Time taken: O(|a| + [M/w]).
 *
 * @param a A substring
 * @param B2 A WordVector to AND with
 * @return True if a was found in suffix tree and B2 was updated, otherwise false
 */
bool MultiEDSM::occVector(const string & a, WordVector & B2)
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
        unsigned int nodeId = this->STp.id(explicitNode);
        if (this->OVMemU8.find(nodeId) != this->OVMemU8.end())
        {
            this->WordVectorAND_IP(B2, this->OVMemU8[nodeId]);
            return true;
        }
        else if (this->OVMem8.find(nodeId) != this->OVMem8.end())
        {
            unsigned int i, k, sn;
            if (this->STp.is_leaf(explicitNode))
            {
                i = this->OVMem8[nodeId].start;
                sn = this->STp.csa[i];
                if (
                    (sn > 0 && ((sn + 1) < (this->R - 1))) && \
                    (this->STpIdx2BVIdx[sn-1] != SEPARATOR_DIGIT) && \
                    (this->STpIdx2BVIdx[sn]   != SEPARATOR_DIGIT) && \
                    (this->STpIdx2BVIdx[sn+1] != SEPARATOR_DIGIT)
                ) {
                    k = this->STpIdx2BVIdx[sn] - 1;  //calculate bit position -- k is index in bitvector
                    unsigned int wordIdx = (unsigned int) ((double)k / (double)BITSINWORD);
                    unsigned int bitShifts = k % BITSINWORD;
                    WORD h = 1ul << bitShifts;
                    if (B2[wordIdx] & h)
                    {
                        for (i = 0; i < B2.size(); i++) {
                            B2[i] = 0ul;
                        }
                        this->WordVectorSet1At(B2, k);
                        return true;
                    }
                }
                return false;
            }
            else
            {
                WordVector v;
                bool oneSet = false;
                for (i = this->OVMem8[nodeId].start; i <= this->OVMem8[nodeId].end; i++)
                {
                    sn = this->STp.csa[i];
                    if (
                        (sn > 0 && ((sn + 1) < (this->R - 1))) && \
                        (this->STpIdx2BVIdx[sn-1] != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn]   != SEPARATOR_DIGIT) && \
                        (this->STpIdx2BVIdx[sn+1] != SEPARATOR_DIGIT)
                    ) {
                        k = this->STpIdx2BVIdx[sn] - 1;  //calculate bit position -- k is index in bitvector
                        this->WordVectorSet1At(v, k);
                        oneSet = true;
                    }
                }
                if (oneSet) {
                    this->WordVectorAND_IP(B2, v);
                    return true;
                }
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
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
 *  - O([M/w]) worst-case time for step 3 bitwise-OR operation and O(min((m + [M/w]), (m([M/w]))) for LeftShift
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
    unsigned int i, m;
    Segment::const_iterator stringI;

    // cout << "Position #" << (this->d + this->D) << endl;

    //search first segment, priming B, then search the next segments normally in the else condition down below
    if (!this->primed)
    {
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
        this->buildBorderPrefixWordVector(S, this->B);

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
        this->buildBorderPrefixWordVector(S, this->B1);

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
                    this->WordVectorOR_IP(this->B1, this->B);
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

                //
                //step 1 - if string is longer than pattern then do a full search on it
                //
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

                //
                //step 2 - if string is a suffix of a previously determined prefix or infix, then report a match
                //
                if (this->B2.size() != this->B.size()) {
                    this->B2 = this->B;
                } else {
                    for (i = 0; i < this->B.size(); i++) {
                        this->B2[i] = this->B[i];
                    }
                }
                if (this->umsa->searchOnState(*stringI, this->B2, 0, this->maxP - 1))
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

                //
                //step 3 - if string is short, then it could be an infix, so determine it is actual infix and mark its position accordingly
                //
                if (m < (this->maxP - 1))
                {
                    this->Np++;
                    this->Nm += m;
                    for (i = 0; i < this->B.size(); i++) {
                        this->B2[i] = this->B[i];
                    }
                    if (this->occVector(*stringI, this->B2))
                    {
                        this->WordVectorLeftShift_IP(this->B2, m);
                        this->WordVectorOR_IP(this->B1, this->B2);
                    }
                }
            }
        }

        //update the state for the next segment search
        for (i = 0; i < this->B1.size(); i++) {
            this->B[i] = this->B1[i];
        }

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
 * The bitwise OR operation performed 'in place' on parameter a. Time taken:
 * O([M/w])
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
        a[i] = 0ul;
    }
}

/**
 * A special left shifting algorithm for use with OccVector result. It ensures a
 * shift operation is not performed if it overlaps into the next pattern in the
 * bitvector. This method uses the Pos2PatId datastructure to identify if a shift
 * illegally crosses into the next pattern, so requires O(M) space and then it
 * saves the shifted bits temporarily into a WordVector, so the final space usage
 * is O(M + [M/w]).
 * Time taken is roughly: O(num_set_bits(x) + [M/w])... or really
 * worst-Worst-WORST case, where there is a match at every single position, which
 * is very-Very-VERY unlikely then O(M + [M/w])!
 *
 * @param x
 * @param m The length of the string being checked
 */
void MultiEDSM::WordVectorSPECIALSHIFT_IP(WordVector & x, unsigned int m)
{
    WORD temp;
    unsigned int j, currPos, updPos, currPattId, updPattId, currWordIdx, updWordIdx, n = x.size();
    WordVector y(n, 0ul);
    int i = n - 1;
    while (i > -1)
    {
        temp = x[i];
        while (temp)
        {
            j = ffs(temp) - 1;
            currPos = i * BITSINWORD + j;
            currPattId = this->Pos2PatId[currPos];
            updPos = currPos + m;
            updPattId = this->Pos2PatId[updPos];
            if (updPos < this->M && currPattId == updPattId)
            {
                updWordIdx = (unsigned int) ((double)updPos / (double)BITSINWORD);
                y[updWordIdx] = y[updWordIdx] | (1ul << (updPos % BITSINWORD));
            }
            temp = temp ^ (1ul << j);
        }
        i--;
    }
    for (j = 0; j < n; j++) {
        x[j] = y[j];
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
    unsigned int i = 0, p = 0, s = x.size();
    unsigned int specialShiftTime = 0, simpleShiftTime = m * s;

    while (i < s && specialShiftTime <= simpleShiftTime)
    {
        p += popcount(x[i]);
        specialShiftTime = p + s;
        i++;
    }

    if (simpleShiftTime <= specialShiftTime) {
        this->WordVectorSIMPLESHIFT_IP(x, m);
    } else {
        this->WordVectorSPECIALSHIFT_IP(x, m);
    }
}

/**
 * Sets the bit at position pos to 1 in WordVector x. 0-based bit index. If x is
 * too small it will be expanded to accomodate the new bit.
 *
 * @param x WordVector
 * @param pos The position to place the 1 at. 0-based bit index.
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
