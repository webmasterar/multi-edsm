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

    cout << "Got here 0.1" << endl;

    //add patterns to bitvector and build suffix tree string
    string p = "";
    char sep = 0;
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
        while (this->alphabet.find(sep) != string::npos);

        //update STpIdx2BVIdx datastructure with correct index of STp match for bitvector
        h = this->patterns[i].length();
        for (j = 0; j < (h + 1); j++) {
            this->STpIdx2BVIdx.push_back(this->M + j - i);
        }
        this->M += h;

        //find min/max pattern length
        if (h < this->minP) {
            this->minP = h;
        }
        if (h > this->maxP) {
            this->maxP = h;
        }
    }

    cout << "Got here 0.2" << endl;

    //construct the suffix tree of patterns with seperators
    construct_im(this->STp, p.c_str(), sizeof(char));

    cout << "Got here 0.3" << endl;

    //stop timer
    this->duration += clock() - start;
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
    bool isDeterminateSegment = false;
    unsigned int matchIdx, posIdx, pattId;

    //temp state variables, segment string iterator
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

        cout << "Got here 1" << endl;

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
                                cout << "Got here 2" << endl;
                                posIdx = match.first;
                                matchIdx = this->pos + posIdx + 1;
                                pattId = match.second;
                                this->report(matchIdx, posIdx, pattId);
                                //@TODO it might be that the same posIdx matches multiple patterns so just report the first one?
                            }
                        }
                        else
                        {
                            if (this->reportOnce)
                            {
                                const pair<int,int> match = *this->umsa->getMatches().begin();
                                posIdx = match.first;
                                matchIdx = this->pos + 1;
                                pattId = match.second;
                                this->report(matchIdx, posIdx, pattId);
                                break;
                            }
                            else
                            {
                                for (const pair<int,int> & match : this->umsa->getMatches()) {
                                    posIdx = match.first;
                                    matchIdx = this->pos + 1;
                                    pattId = match.second;
                                    this->report(matchIdx, posIdx, pattId);
                                }
                            }
                        }
                        matchFound = true;
                    }
                }
            }
        }

        //build equivalent of BorderPrefixTable
        this->B = this->buildBorderPrefixWordVector(S);

        this->primed = true;
    }

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
 * @param a
 * @param b
 * @return (a & b)
 */
WordVector MultiEDSM::WordVectorAND(const WordVector & a, const WordVector & b)
{
    WordVector c;
    unsigned int i, m = min(a.size(), b.size()), n = max(a.size(), b.size());
    for (i = 0; i < m; i++) {
        c.push_back(a[i] & b[i]);
    }
    for (; i < n; i++) {
        c.push_back(0ul);
    }
    return c;
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
    WordVector C;
    unsigned int m, i = 0;
    Segment::const_iterator stringI;

    for (stringI = S.begin(); stringI != S.end(); ++stringI)
    {
        if (*stringI != EPSILON)
        {
            m = (*stringI).length();
            if (m >= this->maxP) {
                this->umsa->search((*stringI).substr(m - this->maxP));
            } else {
                this->umsa->search(*stringI);
            }
            if (i == 0) {
                C = this->umsa->getLastSearchState();
            } else {
                C = this->WordVectorOR(C, this->umsa->getLastSearchState());
            }
            i++;
        }
    }

    return C;
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
