//License MIT 2017 Ahmad Retha

#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <map>

#ifndef __UMSA__
#define __UMSA__

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

class UnrestrictedMultiShiftAnd
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
     * @var matches A map of ending positions where a match was found and the id
     * of the pattern discovered, <index_in_t, pattern_id>
     */
    std::map<int,int> matches;

    /**
     * @var positions A map of ending positions and the id of the pattern that is found there
     */
    std::map<int,int> positions;

public:
    UnrestrictedMultiShiftAnd(const std::string & alphabet);
    void addPattern(const std::string & pattern);
    bool search(const std::string & text);
    bool search(const std::string & text, std::vector<WORD> & startingSearchState);
    std::map<int,int> getMatches() const;
    std::vector<WORD> getLastSearchState() const;
    void clearMatches();
    unsigned int getNumberOfPatterns() const;
    unsigned int getTotalPatternLength() const;
};

#endif
