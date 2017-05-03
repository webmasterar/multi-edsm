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
#include "UnrestrictedMultiShiftAnd.hpp"
#include "MultiEDSM.hpp"

using namespace sdsl;
using namespace std;

/**
 * @constructor
 * @param alphabet Alphabet used by the patterns
 * @param patterns The patterns to look for
 */
MultiEDSM::MultiEDSM(const std::string & alphabet, const std::vector<std::string> & patterns)
{
    this->alphabet = alphabet;
    this->patterns = patterns;
    this->minP = UINT_MAX;
    this->maxP = 0;
    this->f = 0;
    this->F = 0;
    this->d = 0;
    this->D = 0;
    this->M = 0;
    this->Np = 0;
    this->Nm = 0;
    this->pos = 0;
    this->primed = false;
    this->duration = 0;
    this->reportOnce = true;
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
 * Proprocess the patterns - initialize the bitvector and the suffix tree tools
 */
void MultiEDSM::preprocessPatterns()
{
    //start timer
    clock_t start = clock();

    //initialize Shift-And searching tool
    this->umsa = new UnrestrictedMultiShiftAnd(this->alphabet);
    this->umsa->reportPatternIds(this->reportPatterns);

    //add patterns to bitvector and build suffix tree string
    string p = "";
    char sep = 1;
    unsigned int h, i, j;
    for (i = 0; i < this->patterns.size(); i++)
    {
        //create bitvector of pattern - for Shift-And
        this->umsa->addPattern(this->patterns[i]);

        //add pattern to p and append with unique seperator - for suffix tree
        p += this->patterns[i] + sep;
        do {
            sep = (sep + 1) % (CHAR_MAX + 1);
        }
        while ((this->alphabet.find(sep) != string::npos) || (sep == 0));

        //update STpIdx2BVIdx datastructure with correct index of STp match for bitvector
        h = this->patterns[i].length();
        for (j = 0; j < h; j++) {
            this->STpIdx2BVIdx.push_back(this->M + j - i);
        }
        this->STpIdx2BVIdx.push_back(SEPARATOR_DIGIT);
        this->M += h + 1;

        //find min/max pattern length
        if (h < this->minP) {
            this->minP = h;
        }
        if (h > this->maxP) {
            this->maxP = h;
        }
    }
    p.pop_back();

    //construct the suffix tree of patterns with seperators
    construct_im(this->STp, p.c_str(), sizeof(char));

    //stop timer
    this->duration += clock() - start;
}

/**
 * The Shift-And algorithm performs the equivalent of creating a BorderPrefixTable
 * and this function just returns a WordVector with the prefixes of patterns
 * identified in the suffixes of all strings in a segment.
 *
 * @param S A segment
 * @return A WordVector with the prefixes marked
 */
WordVector MultiEDSM::buildBorderPrefixWordVector(const Segment & S)
{
    WordVector c;
    unsigned int m, i = 0;
    Segment::const_iterator stringI;

    for (stringI = S.begin(); stringI != S.end(); ++stringI)
    {
        if (*stringI != EPSILON)
        {
            m = (*stringI).length();
            if (m >= this->maxP) {
                this->umsa->search((*stringI).substr(m - this->maxP + 1));
            } else {
                this->umsa->search(*stringI);
            }
            if (i == 0) {
                c = this->umsa->getLastSearchState();
            } else {
                c = this->WordVectorOR(c, this->umsa->getLastSearchState());
            }
            i++;
        }
    }

    return c;
}

/**
 * occVector function for demarcating valid infix starting positions. @TODO rewrite
 *
 * @param a @TODO
 * @return @TODO
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

    WordVector v(1, 0ul);

    //if a is present in P
    if (j == a.length()) {
        this->recFindAllChildNodes(explicitNode, v, j);
    }

    return v;
}

/**
* Recursively finds leaves in the tree from node u and updates v
*
* @param u The starting node
* @param v The bitvector where leaves string lengths are encoded to
* @param j The length of the match to move bits along correctly
*/
void MultiEDSM::recFindAllChildNodes(const cst_node_t & u, WordVector & v, const unsigned int j)
{
    if (this->STp.is_leaf(u))
    {
        //sn(u) gets the suffix index from the leaf node u, while STpIdx2BVIdx converts it to the correct index
        int h = (int)this->STp.sn(u); //h is index in suffix tree string
        int i = (int)this->STpIdx2BVIdx[h]; //i is index in bitvector

        //check that it is actually an infix (does not start at position 0 or end with a separator)
        if ((h - 1) >= 0 && this->STpIdx2BVIdx[h - 1] != SEPARATOR_DIGIT && (h + j) < this->STpIdx2BVIdx.size() && this->STpIdx2BVIdx[h + j] != SEPARATOR_DIGIT)
        {
            int wordIdx = (int) ((float)i / (float)BITSINWORD);

            //check if it's a valid infix starting position
            if (this->B[wordIdx] & (1 << ((i % BITSINWORD) - 1)))
            {
                //update position
                i = i + j;
                wordIdx = (int) ((float)i / (float)BITSINWORD);

                //check the word vector has enough words in it
                if (wordIdx >= (int)v.size()) {
                    int k;
                    for (k = v.size(); k <= wordIdx; k++) {
                        v.push_back(0ul);
                    }
                }

                //write position to v
                v[wordIdx] = v[wordIdx] | (1 << ((i % BITSINWORD) - 1));
            }
        }
    }
    else
    {
        for (const auto & child : this->STp.children(u)) {
            this->recFindAllChildNodes(child, v, j);
        }
    }
}

/**
 * Report match found. The results can be obtained by calling MultiEDSM::getMatches()
 *
 * @param matchIdx The index of the ending position where the match was found
 * @param posIdx The position inside a segment where the match was found
 * @param pattId The pattern id
 */
void MultiEDSM::report(const unsigned int matchIdx, const unsigned int posIdx, const unsigned int pattId)
{
    cout << "Match found. matchIdx " << matchIdx << ", posIdx " << posIdx << ", pattId " << pattId << endl;
}

/**
 * Search for the patterns in a segment
 */
bool MultiEDSM::searchNextSegment(const Segment & S)
{
    //start timer
    clock_t start = clock();

    //match variables
    bool matchFound = false;
    bool suffixMatchFound = false;
    bool isDeterminateSegment = false;
    unsigned int matchIdx, posIdx, pattId;

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
                                //@TODO it might be that the same posIdx matches multiple patterns so just report the first one?
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
            this->pos = this->f;
        } else {
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

                //step 1 - if string is longer than patten then do a full search on it
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
                if (m < this->maxP) {
                    suffixMatchFound = this->umsa->search(*stringI, B2);
                } else {
                    suffixMatchFound = this->umsa->search(*stringI, B2, m - this->maxP);
                }
                if (suffixMatchFound) {
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
                    matchFound = true;
                }

                //step 3 - if string is short, then it could be an infix, so determine it is actual infix and mark its position accordingly
                if (m < (this->maxP - 1))
                {
                    this->Np++;
                    this->Nm += m;
                    B2 = this->occVector(*stringI);
                    B1 = this->WordVectorOR(B1, B2);
                }
            }
        }

        this->B = B1;

        //update position
        if (isDeterminateSegment) {
            this->pos += S[0].length();
        } else {
            this->pos++;
        }
    }

    //stop timer
    this->duration += clock() - start;

    return matchFound;
}

/**
 * The bitwise OR operation performed on WordVectors.
 *
 * @param a
 * @param b
 * @return (a | b)
 */
WordVector MultiEDSM::WordVectorOR(const WordVector & a, const WordVector & b)
{
    WordVector c;
    unsigned int i, m = min(a.size(), b.size());
    for (i = 0; i < m; i++) {
        c.push_back(a[i] | b[i]);
    }
    return c;
}

/**
 * The bitwise AND operation performed on WordVectors.
 *
 * @deprecated Unused, consider deleting @TODO
 * @param a
 * @param b
 * @return (a & b)
 */
WordVector MultiEDSM::WordVectorAND(const WordVector & a, const WordVector & b)
{
    WordVector c;
    unsigned int i;
    unsigned int m = min(a.size(), b.size());
    unsigned int n = max(a.size(), b.size());
    for (i = 0; i < m; i++) {
        c.push_back(a[i] & b[i]);
    }
    for (; i < n; i++) {
        c.push_back(0ul);
    }
    return c;
}

/**
* Get a list of the positions of matches found.
*
* @return @TODO
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
