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
    this->D.push_back(0ul);
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
    unsigned int i, j, m = pattern.length();

    //if we need more memory to hold the new pattern then create and initialize new words
    unsigned int l = (unsigned int) ceil((double)(this->M + m) / (double)BITSINWORD);
    unsigned int numWordsDiff = l - this->L;
    if (numWordsDiff > 0)
    {
        for (i = 0; i < numWordsDiff; i++)
        {
            this->Sv.push_back(0ul);
            this->Ev.push_back(0ul);
            this->D.push_back(0ul);
        }
        for (i = 0; i < numWordsDiff; i++)
        {
            vector<WORD> v;
            for (j = 0; j < this->alphabet.length(); j++)
            {
                v.push_back(0ul);
            }
            this->Bv.push_back(v);
        }
    }

    //process the pattern for Bitvector Bv
    int charIdx;
    unsigned int currWordIdx = (unsigned int) ((float)this->M / (float)BITSINWORD);
    unsigned int currBitIdx = this->M % BITSINWORD;
    for (i = 0; i < m; i++)
    {
        //keep a record of pattern id at this position - for search results
        this->positions.push_back(this->N);

        //mark the position of characters in the Bv bitvector
        charIdx = (int) this->Sigma[(int)pattern[i]];
        if (charIdx > 0) {
            this->Bv[currWordIdx][charIdx - 1] = this->Bv[currWordIdx][charIdx - 1] | (1ul << currBitIdx);
        }
        currBitIdx = ++currBitIdx % BITSINWORD;
        if (currBitIdx == 0) {
            currWordIdx++;
        }
    }

    //process the pattern to set Start bit in Sv
    currWordIdx = (unsigned int) ((float)this->M / (float)BITSINWORD);
    currBitIdx = this->M % BITSINWORD;
    this->Sv[currWordIdx] = this->Sv[currWordIdx] | (1ul << currBitIdx);

    //process the pattern to set End bit in Ev and update the class variables
    this->L = this->L + numWordsDiff;
    this->M = this->M + m;
    this->N = this->N + 1;
    this->Ev[this->L - 1] = this->Ev[this->L - 1] | (1ul << ((this->M % BITSINWORD) - 1));
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
 * Search a text starting at position i
 *
 * @param text
 * @param pos position in text to start searching from
 * @param len The number of characters to search from position pos
 * @return True if one or more matches found, otherwise False
 */
bool MyUMSA::search(const string & text, unsigned int pos, unsigned int len)
{
    //initialize vector D to have a fresh state for the new search
    this->D.assign(this->L, 0ul);
    return this->search(text, this->D, pos, len);
}

/**
 * Search a text but supply an initial search state - this is useful for searching
 * text provided intermittently (on-line algorithm).
 *
 * @param text
 * @param startingSearchState A vector<WORD> with L elements
 * @return True if one or more matches found, otherwise False
 */
bool MyUMSA::search(const string & text, vector<WORD> & startingSearchState)
{
    return this->search(text, startingSearchState, 0, text.length());
}

/**
 * Search a text but supply an initial search state - this is useful for searching
 * text provided intermittently (online algorithm). For simple searches just use
 * MyUMSA::search()
 *
 * @param text
 * @param startingSearchState A vector<WORD> with L elements
 * @param pos position in text to start searching from
 * @return True if one or matches found, otherwise False
 */
bool MyUMSA::search(const string & text, vector<WORD> & startingSearchState, unsigned int pos)
{
    return this->search(text, startingSearchState, pos, text.length());
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
bool MyUMSA::search(const string & text, vector<WORD> & startingSearchState, unsigned int pos, unsigned int len)
{
    unsigned int i, j, k, n;

    if (this->M == 0) {
        cerr << "No patterns added to search." << endl;
        return false;
    }

    //Make sure vector D has sufficient memory for the search
    this->D = startingSearchState;
    if (this->D.size() < this->L)
    {
        j = this->L - this->D.size();
        for (i = 0; i < j; i++) {
            this->D.push_back(0ul);
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
        carry = 0;
        charIdx = (int) this->Sigma[(int)text[i]] - 1;

        if (charIdx == -1)
        {
            if (!zeroed)
            {
                //character not in patterns so clear D
                for (j = 0; j < this->L; j++)
                {
                    this->D[j] = 0;
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

                carry = (WORD) ((carryMask & temp) != 0);
            }
            zeroed = false;
        }
    }

    return matchFound;
}

/**
 * After doing a search, get the last state of the search
 */
vector<WORD> MyUMSA::getLastSearchState() const
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
vector<WORD> MyUMSA::getEndingStates() const
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
vector<unsigned int> MyUMSA::getPatternPositions() const
{
    return this->positions;
}
