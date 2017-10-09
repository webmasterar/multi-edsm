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
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include "MyUMSA.hpp"

using namespace std;

/**
 * @constructor An alphabet must be supplied of the letters that will be used in the patterns
 *
 * @param alphabet A list of characters used in the text/patterns e.g. ACGT
 * @param reportPatterns Report positions when matches found?
 */
MyUMSA::MyUMSA(const string & alphabet, bool reportPatterns)
{
    this->alphabet = alphabet;
    this->reportPatterns = reportPatterns;
    this->N = 0;
    this->M = 0;
    this->L = 1;
    this->Sv.push_back(0ul);
    this->Ev.push_back(0ul);
    //this->D.push_back(0ul);
    vector<WORD> v;
    unsigned int i;
    for (i = 0; i < this->alphabet.length(); i++)
    {
        this->Sigma[(int)alphabet[i]] = i + 1;
        v.push_back(0ul);
    }
    this->Bv.push_back(v);
}

/**
 * Add a pattern
 *
 * @param pattern
 */
void MyUMSA::addPattern(const string & pattern)
{
    int charIdx;
    unsigned int wordIdx, bitIdx, i, j, m = pattern.length();
    unsigned int lastIdx = m - 1;

    for (i = 0; i < m; i++)
    {
        charIdx = (int) this->Sigma[(int)pattern[i]];
        wordIdx = (unsigned int) ((double)(this->M + i) / (double)BITSINWORD);
        bitIdx = (unsigned int) ((this->M + i) % BITSINWORD);

        //expand memory if required
        if (this->L == wordIdx)
        {
            this->Sv.push_back(0ul);
            this->Ev.push_back(0ul);
            //this->D.push_back(0ul);
            vector<WORD> v;
            for (j = 0; j < this->alphabet.length(); j++)
            {
                v.push_back(0ul);
            }
            this->Bv.push_back(v);
            this->L = this->L + 1;
        }

        //keep a record of pattern id at this position - for search results
        this->positions.push_back(this->N);
        //mark letter position of pattern
        if (charIdx > 0) {
            this->Bv[wordIdx][charIdx - 1] = this->Bv[wordIdx][charIdx - 1] | (1ul << bitIdx);
        }
        //mark start position of pattern
        if (i == 0) {
            this->Sv[wordIdx] = this->Sv[wordIdx] | (1ul << bitIdx);
        }
        //mark end position of pattern
        if (i == lastIdx) {
            this->Ev[wordIdx] = this->Ev[wordIdx] | (1ul << bitIdx);
        }
    }

    this->M = this->M + m;
    this->N = this->N + 1;
}

/**
 * Search a text
 *
 * @param text
 * @return True if one or matches found, otherwise False
 */
bool MyUMSA::search(const string & text)
{
    return this->search(text, 0, text.length());
}

/**
 * Search a text starting at position i
 *
 * @param text
 * @param pos position in text to start searching from
 * @return True if one or more matches found, otherwise False
 */
bool MyUMSA::search(const string & text, unsigned int pos)
{
    return this->search(text, pos, text.length());
}

/**
 * Search a text starting at position pos for len characters
 *
 * @param text
 * @param pos position in text to start searching from
 * @param len The number of characters to search from position pos
 * @return True if one or matches found, otherwise False
 */
bool MyUMSA::search(const string & text, unsigned int pos, unsigned int len)
{
    unsigned int i, j, k, n;

    if (this->M == 0) {
        cerr << "No patterns added to search." << endl;
        return false;
    }

    //initialize vector D to have a fresh state for the new search
    if (this->D.size() != this->L) {
        this->D.assign(this->L, 0ul);
    } else {
        for (i = 0; i < this->L; i++) {
            this->D[i] = 0ul;
        }
    }

    //calculate starting and ending positions of search
    i = min((int)pos, (int)text.length() - 1);
    n = min(i + len, (unsigned int)text.length());

    //init tracking vars
    int charIdx;
    bool zeroed = false;
    bool matchFound = false;
    WORD temp, carry, check;
    WORD carryMask = 1ul << (BITSINWORD - 1);

    //loop through the text
    for (; i < n; i++)
    {
        carry = 0ul;
        charIdx = (int) this->Sigma[(int)text[i]] - 1;

        if (charIdx == -1)
        {
            if (!zeroed)
            {
                //character not in patterns so clear D
                for (j = 0; j < this->L; j++)
                {
                    this->D[j] = 0ul;
                }
                zeroed = true;
            }
        }
        else
        {
            //loop through the words
            for (j = 0; j < this->L; j++)
            {
                temp = this->D[j];

                this->D[j] = (((this->D[j] << 1) | carry) | this->Sv[j]) & this->Bv[j][charIdx];

                //check if any matches found
                if (this->reportPatterns)
                {
                    check = this->D[j] & this->Ev[j];
                    while (check)
                    {
                        matchFound = true;
                        k = ffs(check) - 1;
                        this->matches.push_back(pair<int,int>((int)i, this->positions[j * BITSINWORD + k]));
                        check = check ^ (1ul << k);
                    }
                }

                carry = (WORD) ((carryMask & temp) != 0ul);
            }
            zeroed = false;
        }
    }

    return matchFound;
}

/**
 * Search a text but supply an initial search state - this is useful for searching
 * text provided intermittently (online algorithm). For simple searches just use
 * MyUMSA::search()
 *
 * @param text
 * @param startingSearchState A vector<WORD> with L elements
 * @param pos position in text to start searching from
 * @param len The number of characters to search from position pos
 * @return True if one or matches found, otherwise False
 */
bool MyUMSA::searchOnState(const string & text, vector<WORD> & startingSearchState, unsigned int pos, unsigned int len)
{
    unsigned int i, j, k, n;

    if (this->M == 0) {
        cerr << "No patterns added to search." << endl;
        return false;
    }

    //Make sure vector startingSearchState has sufficient memory for the search
    if (startingSearchState.size() < this->L)
    {
        j = this->L - startingSearchState.size();
        for (i = 0; i < j; i++) {
            startingSearchState.push_back(0ul);
        }
    }

    //calculate starting and ending positions of search
    i = min((int)pos, (int)text.length() - 1);
    n = min(i + len, (unsigned int)text.length());

    //init tracking vars
    int charIdx;
    bool zeroed = false;
    bool matchFound = false;
    WORD temp, carry, check;
    WORD carryMask = 1ul << (BITSINWORD - 1);

    //loop through the text
    for (; i < n; i++)
    {
        carry = 0ul;
        charIdx = (int) this->Sigma[(int)text[i]] - 1;

        if (charIdx == -1)
        {
            if (!zeroed)
            {
                //character not in patterns so clear D
                for (j = 0; j < this->L; j++)
                {
                    startingSearchState[j] = 0ul;
                }
                zeroed = true;
            }
        }
        else
        {
            //loop through the words
            for (j = 0; j < this->L; j++)
            {
                temp = startingSearchState[j];

                startingSearchState[j] = (((startingSearchState[j] << 1) | carry) | this->Sv[j]) & this->Bv[j][charIdx];

                //check if any matches found
                if (this->reportPatterns)
                {
                    check = startingSearchState[j] & this->Ev[j];
                    while (check)
                    {
                        matchFound = true;
                        k = ffs(check) - 1;
                        this->matches.push_back(pair<int,int>((int)i, this->positions[j * BITSINWORD + k]));
                        check = check ^ (1ul << k);
                    }
                }

                carry = (WORD) ((carryMask & temp) != 0ul);
            }
            zeroed = false;
        }
    }

    return matchFound;
}

/**
 * After doing a search, get the last state of the search
 *
 * @return The last state of the search
 */
const vector<WORD> & MyUMSA::getLastSearchState() const
{
    return this->D;
}

/**
 * Clear matches
 */
void MyUMSA::clearMatches()
{
    this->matches.clear();
}

/**
 * Get the number of patterns already added
 */
unsigned int MyUMSA::getNumberOfPatterns() const
{
    return this->N;
}

/**
 * Get the total length of the patterns added
 */
unsigned int MyUMSA::getTotalPatternLength() const
{
    return this->M;
}

/**
 * Get the length of the state vector (how many computer words are used)
 */
unsigned int MyUMSA::getStateVectorLength() const
{
    return this->L;
}

/**
 * Get a multimap of all the matches found - <index_in_t, pattern_id>
 */
vector<pair<int,int>> MyUMSA::getMatches() const
{
    return this->matches;
}

/**
 * This function simply returns a WordVector containing the ending states of
 * the Shift-And patterns, where the 1s are indicators of matching positions
 * for when a pattern successfully matches in a searched text
 *
 * @return A WordVector marking the ending positions/states of the pattern set
 */
const vector<WORD> & MyUMSA::getEndingStates() const
{
    return this->Ev;
}

/**
 * Enable reporting of matches
 */
void MyUMSA::enableReporting()
{
    this->reportPatterns = true;
}

/**
 * Disable reporting of matches
 */
void MyUMSA::disableReporting()
{
    this->reportPatterns = false;
}

/**
 * Returns a list where for any given position it contains the id of the pattern
 *
 * @return vector containing pattern id at every position in the bitvector
 */
const vector<unsigned int> & MyUMSA::getPatternPositions() const
{
    return this->positions;
}
