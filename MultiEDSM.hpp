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

#ifndef __MULTI_EDSM__
#define __MULTI_EDSM__

#include <cstdlib>
#include <string>
#include <vector>
#include <memory>
#include <sdsl/util.hpp>
#include <sdsl/suffix_trees.hpp>
#include "UnrestrictedMultiShiftAnd.hpp"

#define EPSILON "E"
#define BUFFERSIZE 1000000
#define SEPARATOR_DIGIT -1

typedef sdsl::cst_sct3<> cst_t;
typedef cst_t::node_type cst_node_t;
typedef cst_t::size_type cst_size_t;
typedef cst_t::char_type cst_char_t;
typedef std::vector<std::tuple<unsigned int, unsigned int>> ResultSet;
typedef std::vector<std::string> Segment;
typedef std::vector<WORD> WordVector;
typedef std::vector<Segment> GenIndSeq;

class MultiEDSM
{
private:

    void recFindAllChildNodes(const cst_node_t & u, WordVector & v, const unsigned int j);

    void preprocessPatterns();

    WordVector WordVectorOR(const WordVector & a, const WordVector & b);

    WordVector WordVectorAND(const WordVector & a, const WordVector & b);

protected:
    /**
     * @var matches All matches found as tuple of <ending_position, pattern_id>
     */
    ResultSet matches;

    /**
     * @var alphabet The alphabet of the patterns
     */
    std::string alphabet;

    /**
     * @var patterns The determinate patterns to find in T
     */
    std::vector<std::string> patterns;

    /**
     * @var STpIdx2BVIdx A tool to return the correct index of the bitvector from the suffix tree position - uses O(M) space, O(1) time
     */
    std::vector<int> STpIdx2BVIdx;

    /**
     * @var M The total length of the patterns
     */
    unsigned int M;

    /**
     * @var minP The length of the smallest pattern in patterns
     */
    unsigned int minP;

    /**
     * @var maxP The length of the longest pattern in patterns
     */
    unsigned int maxP;

    /**
     * @var pos The current position in the input passed in. Degenerate
     * positions/segments count as 1 position, whereas determinate segments are
     * counted as N positions (as many characters as there are in the segment).
     */
    unsigned int pos;

    /**
     * @var duration The amount of time spent by EDSM-BV
     */
    double duration;

    /**
     * @var primed Has the algorithm been primed with an initial segment to search?
     */
    bool primed;

    /**
     * @var reportOnce Only report once for a position in indeterminate segments -- default: True
     */
    bool reportOnce;

    /**
     * @var reportPatterns Report id of pattern found -- default: True
     */
    bool reportPatterns;

    /**
     * @var STp The suffix tree (array) of P
     */
    cst_t STp;

    /**
     * @var B The current state of the search
     */
    WordVector B;

    /**
     * @var umsa UnrestrictedMultiShiftAnd object used for searching - stores the bitvector of the patterns
     */
    UnrestrictedMultiShiftAnd * umsa;

    /**
     * @var f The total length of determinate segments searched
     */
    unsigned int f;

    /**
     * @var F the total length of degenerate segments searched
     */
    unsigned int F;

    /**
     * @var d The number of determinate segments searched so far
     */
    unsigned int d;

    /**
     * @var D The number of degenerate segments searched so far
     */
    unsigned int D;

    /**
     * @var Np The total number of strings analyzed with length < minP
     */
    unsigned int Np;

    /**
     * @var Nm The total length of strings analyzed with length < minP
     */
    unsigned int Nm;

    void report(const unsigned int matchIdx, const unsigned int posIdx, const unsigned int pattId);

    WordVector buildBorderPrefixWordVector(const Segment & S);

    WordVector occVector(const std::string & a);

public:

    MultiEDSM(const std::string & alphabet, const std::vector<std::string> & patterns);

    ~MultiEDSM();

    bool searchNextSegment(const Segment & S);

    ResultSet getMatches() const;

    void reportOncePerPosition(bool yesorno = true);

    void reportPatternIds(bool yesorno = true);

    void clearMatches();

    double getDuration() const;

    unsigned int getd() const;

    unsigned int getD() const;

    unsigned int getf() const;

    unsigned int getF() const;

    unsigned int getNp() const;

    unsigned int getNm() const;

};

#endif
