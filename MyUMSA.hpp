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

#ifndef __MYUMSA__
#define __MYUMSA__

#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <map>

#define WORD unsigned long int
#define WORDSIZE sizeof(WORD)
#define BITSINWORD (WORDSIZE * 8)
#if INTPTR_MAX == INT32_MAX
    #define ffs(x) __builtin_ffs((x))
#elif INTPTR_MAX == INT64_MAX
    #define ffs(x) __builtin_ffsl((x))
#else
    #error "Unsupported architecture - neither 64 or 32 bit!"
#endif

class MyUMSA
{
protected:
    /**
     * @var L The number of computer words needed to store the patterns
     */
    unsigned int L;
    /**
     * @var M The total length of the patterns
     */
    unsigned int M;
    /**
     * @var n The number of patterns added
     */
    unsigned int N;
    /**
     * @var Sv The bitvector holding the initial starting positions (states) of patterns
     */
    std::vector<WORD> Sv;
    /**
     * @var Ev The bitvector holding the ending positions (states) of patterns
     */
    std::vector<WORD> Ev;
    /**
     * @var Bv The bitvector that holds the positions of a character for all
     * characters. The first vector is the WORD position and the second is a
     * vector of WORDS encoding characters and their positions.
     */
    std::vector<std::vector<WORD>> Bv;
    /**
     * @var D Delta, the vector holding the current state of the search
     */
    std::vector<WORD> D;
    /**
     * @var Sigma The indexes of letters using the ascii table of characters
     */
    unsigned char Sigma[256] = {0};
    /**
     * @var alphabet The alphabet to use
     */
    std::string alphabet;
    /**
     * @var matches A multimap of ending positions where a match was found and
     * the id of the pattern discovered, <index_in_t, pattern_id>
     */
    std::vector<std::pair<int,int>> matches;
    /**
     * @var positions A vector holding the id of the pattern found at this position
     * in the bitvector. Although we can use a map, it would take O(N log N) time
     * to construct it and O(log N) time for each look up. This vector will take
     * O(M) space and O(1) time to lookup a match so that's why we use it.
     */
    std::vector<unsigned int> positions;
    /**
     * @var reportPatterns
     */
    bool reportPatterns;

public:
    MyUMSA(const std::string & alphabet, bool reportPatterns = true);
    bool search(const std::string & text);
    bool search(const std::string & text, unsigned int pos);
    bool search(const std::string & text, unsigned int pos, unsigned int len);
    bool search(const std::string & text, std::vector<WORD> & startingSearchState);
    bool search(const std::string & text, std::vector<WORD> & startingSearchState, unsigned int pos);
    bool search(const std::string & text, std::vector<WORD> & startingSearchState, unsigned int pos, unsigned int len);
    std::vector<std::pair<int,int>> getMatches() const;
    std::vector<WORD> getLastSearchState() const;
    std::vector<WORD> getEndingStates() const;
    unsigned int getNumberOfPatterns() const;
    unsigned int getTotalPatternLength() const;
    void addPattern(const std::string & pattern);
    void clearMatches();
    void disableReporting();
    void enableReporting();
};

#endif
