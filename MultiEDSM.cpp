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
#include <ctime>
#include <divsufsort64.h>
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
    this->patterns = patterns;
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
    this->preprocessPatterns();
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
 */
void MultiEDSM::preprocessPatterns()
{
    //start timer
    clock_t start = clock();

    //initialize Shift-And searching tool
    this->umsa = new MyUMSA(this->alphabet, this->reportPatterns);

    //add patterns to bitvector and build suffix tree string
    string p = "";
    char sep = '#';
    unsigned int h, i, j;
    for (i = 0; i < this->patterns.size(); i++)
    {
        h = this->patterns[i].length();
        if (h == 0) {
            cerr << "Error: zero-length pattern given!" << endl;
            throw 1;
        }

        //create bitvector of pattern - for Shift-And
        this->umsa->addPattern(this->patterns[i]);

        //add pattern to p and append with seperator for building suffix tree of patterns
        p += this->patterns[i];
        p += sep;

        //update STpIdx2BVIdx and Pos2PatId datastructures with correct index of
        //STp match for bitvector and pattern id based on bit position in bitvector
        this->M += h;
        for (j = 0; j < h; j++)
        {
            if (this->alphabet.find(patterns[i][j]) == string::npos)
            {
                cerr << "Error: Invalid character found in pattern: " << (int) patterns[i][j] << endl;
                throw 1;
            }
            this->STpIdx2BVIdx.push_back(this->R + j - i);
            this->Pos2PatId.push_back(i);
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
    construct_im(this->STp, p.c_str(), sizeof(char));

    //construct occVector datastructure
    this->constructOV();

    //stop timer
    this->duration += clock() - start;
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
    for (i = 0; i < numNodes; i++) {
        WordVector v(1, 0ul);
        this->OVMem.push_back(v);
    }
    this->recAssignOVMem(this->STp.root());
}

/**
 * Recursively traverses the suffix tree STp and denotes starting positions of
 * substrings, encoding them into WordVectors, and combining them as it moves
 * back up the tree. Space efficient because it only uses as many computer words
 * as is necessary.
 */
WordVector MultiEDSM::recAssignOVMem(const cst_node_t & u)
{
    int id = this->STp.id(u);

    if (this->STp.is_leaf(u))
    {
        int h = (int)this->STp.sn(u); //h is index in suffix tree string
        int i = (int)this->STpIdx2BVIdx[h]; //i is index in bitvector
        int j = this->OVMem[id].size(); //current size of the OVMem[id] WordVector

        //make sure it is not first letter or last letter in a pattern
        if (h > 0 && this->STpIdx2BVIdx[h - 1] != SEPARATOR_DIGIT && \
           (h + 1) < ((int)this->R - 1) && this->STpIdx2BVIdx[h + 1] != SEPARATOR_DIGIT)
        {
            int wordIdx = (int) ((float)(i - 1) / (float)BITSINWORD);
            int bitShifts = (i % BITSINWORD) - 1;

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
            this->OVMem[id] = this->WordVectorOR(this->OVMem[id], this->recAssignOVMem(child));
        }
    }

    return this->OVMem[id];
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
        if (*stringI != EPSILON)
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
                c = this->WordVectorOR(c, this->umsa->getLastSearchState());
            }
            i++;
        }
    }

    this->umsa->enableReporting();

    return c;
}

/**
 * occVector tool returns the starting positions of the given substring a if it
 * is present in the pattern encoded as a bitvector. Time taken: O(|a|).
 *
 * @param a A substring
 * @return Starting positions of substring a encoded into a bitvector
 */
WordVector MultiEDSM::occVector(const string & a)
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

    WordVector v;

    //if a is present in P
    if (j == a.length()) {
        int id = this->STp.id(explicitNode);
        v = this->OVMem[id];
    }

    return v;
}

/**
 * Report match found. The results can be obtained by calling MultiEDSM::getMatches()
 *
 * @param matchIdx The index of the ending position in the EDT where the match was found
 * @param posIdx The position inside a segment where the match was found
 * @param pattId The id of the pattern found
 */
void MultiEDSM::report(const unsigned int matchIdx, const unsigned int posIdx, const int pattId)
{
    cout << "Match found. matchIdx " << matchIdx << ", posIdx " << posIdx << ", pattId " << pattId << endl;
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
    unsigned int m;
    WordVector B1, B2;
    Segment::const_iterator stringI;

    //search first segment, priming B, then search the next segments normally in the else condition
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
                if (*stringI != EPSILON) {
                    this->F += (*stringI).length();
                }
            }
            isDeterminateSegment = false;
        }

        //search for the patterns
        for (stringI = S.begin(); stringI != S.end(); ++stringI)
        {
            if (*stringI != EPSILON)
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
                if (*stringI == EPSILON) {
                    B1 = this->WordVectorOR(B1, this->B);
                } else {
                    this->F += (*stringI).length();
                }
            }
            isDeterminateSegment = false;
        }

        //loop through the strings in the segment
        for (stringI = S.begin(); stringI != S.end(); ++stringI)
        {
            if (*stringI != EPSILON)
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
                    B2 = WordVectorAND(this->B, this->occVector(*stringI));
                    B2 = WordVectorLeftShift(B2, m);
                    B1 = this->WordVectorOR(B1, B2);
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
 * The bitwise OR operation performed on WordVectors. Time taken: O([max_words_used(a,b)/w])
 *
 * @param a
 * @param b
 * @return (a | b)
 */
WordVector MultiEDSM::WordVectorOR(const WordVector & a, const WordVector & b)
{
    WordVector c;
    unsigned int i, m = min(a.size(), b.size()), n = max(a.size(), b.size());
    for (i = 0; i < m; i++) {
        c.push_back(a[i] | b[i]);
    }
    if (a.size() > b.size())
    {
        for (; i < n; i++) {
            c.push_back(a[i]);
        }
    }
    else
    {
        for (; i < n; i++) {
            c.push_back(b[i]);
        }
    }
    return c;
}

/**
 * The bitwise AND operation performed on WordVectors. Time taken: O([min_words_used(a,b)/w])
 *
 * @param a
 * @param b
 * @return (a & b)
 */
WordVector MultiEDSM::WordVectorAND(const WordVector & a, const WordVector & b)
{
    WordVector c;
    unsigned int i;
    unsigned int m = min(a.size(), b.size());
    for (i = 0; i < m; i++) {
        c.push_back(a[i] & b[i]);
    }
    return c;
}

/**
 * A special left shifting algorithm for use with OccVector result. It ensures a
 * shift operation is not performed if it overlaps into the next pattern in the
 * bitvector. This method uses the Pos2PatId datastructure to identify if a shift
 * illegally crosses into the next pattern, so takes O(M) space.
 * Time taken on average I predict: O(num_set_bits(x) + [M/w])... or really
 * worst-Worst-WORST case, where there is a match at every single position, which
 * is very-Very-VERY unlikely then O(M + [M/w])!
 *
 * @param x
 * @param m The length of the string being checked
 * @return A WordVector containing only the states of valid infixes
 */
WordVector MultiEDSM::WordVectorSPECIALSHIFT(const WordVector & x, unsigned int m)
{
    WordVector c;

    WORD temp;
    unsigned int i, j, currPos, updPos, currPattId, updPattId, updWordIdx;
    for (i = 0; i < x.size(); i++)
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
                updWordIdx = (unsigned int) ((float)updPos / (float)BITSINWORD);
                while (updWordIdx >= c.size()) {
                    c.push_back(0ul);
                }
                c[updWordIdx] = c[updWordIdx] | (1ul << (updPos % BITSINWORD));
            }
            temp = temp ^ (1ul << j);
        }
    }

    return c;
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
 * @return A WordVector containing only the states of valid infixes
 */
WordVector MultiEDSM::WordVectorSIMPLESHIFT(const WordVector & x, unsigned int m)
{
    unsigned int s = x.size();
    if (s == 0) {
        return x;
    }
    WordVector c = x;
    WordVector ends = this->umsa->getEndingStates();

    WORD carryMask = 1ul << (BITSINWORD - 1);
    WORD temp, carry;
    unsigned int h, i;
    for (h = 0; h < m; h++)
    {
        carry = 0;
        for (i = 0; i < s; i++)
        {
            temp = c[i];
            c[i] = (c[i] << 1) | carry; //left shift and apply the carried 1
            c[i] = c[i] ^ (ends[i] & c[i]); //if 1 passes pattern boundard (end), then remove it because it is an illegal shift
            carry = (WORD)((carryMask & temp) != 0); //work out if we need to carry a 1 to the next word
            if (carry && (i == (s - 1))) { //identify if we need to expand c
                c.push_back(0ul);
                s++;
            }
        }
    }

    return c;
}

/**
 * Performs the left-shift operation for the WordVector returned from running the
 * OccVector tool. This operation is potentially expensive so there are two
 * different ways to do it:
 *  - MultiEDSM::WordVectorSPECIALSHIFT() has a worst case time of O(num_set_bits(x) + [M/w])
 *    or O(M + [M/w]) if there is a huge number of set bits in x.
 *  - MultiEDSM::WordVectorSIMPLESHIFT() has a worst case time of O(m[M/w])
 *    because it does m shift operations over [M/w] words.
 * So, to get the best performance, we do a triage and identify which of the
 * methods does it in the least time. The triage check itself takes O([M/w]) time.
 *
 * @param x
 * @param m The length of the string being checked
 * @return A WordVector containing only the states of valid infixes
 */
WordVector MultiEDSM::WordVectorLeftShift(const WordVector & x, unsigned int m)
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
        return this->WordVectorSIMPLESHIFT(x, m);
    }
    return this->WordVectorSPECIALSHIFT(x, m);
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
