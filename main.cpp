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

#include <cstdlib>
#include <iostream>
#include <climits>
#include <getopt.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <Variant.h>
#include "MultiEDSM.hpp"

using namespace std;
using namespace vcflib;

struct VarItem
{
    unsigned int pos;
    unsigned int skip;
    Segment seg;
};

typedef vector<struct VarItem> VarItemArray;

char BUFF[BUFFERSIZE];
int BUFFLIMIT = 0;
int POS = 0;

/*
* Buffered file reading char by char from a plain DNA FASTA file
*
* @param f opened file handle
*/
char getNextCharFAS(ifstream & f)
{
    if (BUFFLIMIT == 0) {
        f.read(BUFF, BUFFERSIZE);
        BUFFLIMIT = f.gcount();
        POS = 0;
        if (BUFFLIMIT == 0) {
            return '\0';
        }
    }
    char c = BUFF[POS];
    if (++POS == BUFFLIMIT) {
        BUFFLIMIT = 0;
    }

    if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') {
        return c;
    }
    return getNextCharFAS(f);
}

/*
* Buffered file reading char by char from an EDS file
*
* @param f opened file handle
*/
char getNextCharEDS(ifstream & f)
{
    if (BUFFLIMIT == 0) {
        f.read(BUFF, BUFFERSIZE);
        BUFFLIMIT = f.gcount();
        POS = 0;
        if (BUFFLIMIT == 0) {
            return '\0';
        }
    }
    char c = BUFF[POS];
    if (++POS == BUFFLIMIT) {
        BUFFLIMIT = 0;
    }

    if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' || c == '{' || c == '}' || c == ',' || c == EPSILON[0]) {
        return c;
    }
    return getNextCharEDS(f);
}

/*
* Loops through all the VCF records and generates a master list of all the
* variants in the VCF file. It combines duplicates and nested variants too.
*
* @param vcfName The vcf file to open
* @param variantItems A reference to the array which will be populated with the
* variants (VarItems).
* @return False on error such as VCF file missing
*/
bool populateVarItemArray(string vcfName, VarItemArray & variantItems)
{
    VariantCallFile vf;
    vf.open(vcfName);
    if (!vf.is_open()) {
        cerr << "Error: Failed to open variants file!" << endl;
        return false;
    }
    Variant v(vf);
    unsigned int prevPos = 0, prevRefLen = 0;
    unsigned int currPos, currRefLen;

    while (vf.getNextVariant(v))
    {
        currPos = v.position;
        currRefLen = v.ref.length();
        //handle duplicates
        if (currPos == prevPos)
        {
            //check duplicate doesn't have a longer ref -- if it does merge new ref into segment
            if (currRefLen > prevRefLen && !(v.ref[0] == '<' || v.ref[0] == '.'))
            {
                Segment updSeg;
                for (const string & t : variantItems.back().seg)
                {
                    string r = t;
                    if (r.length() < currRefLen) {
                        r += v.ref.substr(prevRefLen);
                    }
                    updSeg.push_back(r);
                }
                updSeg[0] = v.ref;
                variantItems.back().seg = updSeg;
                variantItems.back().skip = currRefLen;
            }
            //add new alts
            for (const string & t : v.alt)
            {
                if (!(v.ref[0] == '<' || v.ref[0] == '.')) {
                    variantItems.back().seg.push_back(t);
                }
            }
        }
        //handle nested variants
        else if (currPos < (prevPos + prevRefLen))
        {
            Segment lastSeg;
            for (const string & s : variantItems.back().seg)
            {
                for (const string & t : v.alt)
                {
                    if ((prevPos + s.length()) > currPos && !(t[0] == '<' || t[0] == '.')) {
                        string r;
                        r = s.substr(0, currPos - prevPos);
                        r += t;
                        r += s.substr(currPos - prevPos + 1);
                        lastSeg.push_back(r);
                    }
                }
            }
            for (const string & r : lastSeg) {
                variantItems.back().seg.push_back(r);
            }
        }
        //handle regular variants
        else
        {
            Segment currSeg;
            for (const string & t : v.alleles) {
                if (!(t[0] == '<' || t[0] == '.')) {
                    currSeg.push_back(t);
                }
            }
            struct VarItem varItem = {currPos, currRefLen, currSeg};
            variantItems.push_back(varItem);
            prevPos = currPos;
            prevRefLen = currRefLen;
        }
    }
    return true;
}

/**
 * Convert raw string like 3.5g (3.5 gigabytes) to 3584 (megabytes) -- return 0 indicates error
 */
unsigned int getMemLimitMB(string rawString)
{
    int l = rawString.length();
    if (l == 0) {
        return 0;
    }
    char unit = rawString[l - 1];
    if (unit == 'm' || unit == 'M') {
        return (unsigned int) ceil(atof(rawString.substr(0, l - 1).c_str()));
    } else if (unit == 'g' || unit == 'G') {
        return (unsigned int) ceil(atof(rawString.substr(0, l - 1).c_str()) * 1024);
    }
    return 0;
}

/**
 * Search for patterns using the Reference fasta file with variants file
 *
 * @param multiedsm The multiedsm object (reference)
 * @param seqFile The FASTA sequence file name
 * @param varFile The variants file name
 * @return False on error
 */
bool searchVCF(MultiEDSM * multiedsm, string & seqFile, string & varFile)
{
    //open reference file for reading
    ifstream rf(seqFile.c_str(), ios::in);
    if (!rf.good()) {
        cerr << "Error: Failed to open reference file!" << endl;
        return false;
    }

    //create variables for reading through vcf records and looking for duplicates
    VarItemArray variantItems;
    bool parsedVCF = populateVarItemArray(varFile, variantItems);
    if (!parsedVCF) {
        return false;
    }
    VarItemArray::iterator vit = variantItems.begin();
    unsigned int vPos = vit->pos;
    unsigned int vSkip = vit->skip;
    Segment vSeg = vit->seg;

    //initialize fasta file reading and segment creation helper variables
    string tBuff = "";
    tBuff.reserve(BUFFERSIZE);
    char c;
    unsigned int rfIdx = 1, i = 0;
    Segment segment;

    //skip first line of fasta file
    getline(rf, tBuff);
    tBuff = "";

    //go through the reference sequence
    while ((c = getNextCharFAS(rf)) != '\0')
    {
        if (rfIdx != vPos)
        {
            tBuff += c;
            i++;
            if (i >= BUFFERSIZE) {
                segment.clear();
                segment.push_back(tBuff);
                multiedsm->searchNextSegment(segment);
                tBuff = "";
                i = 0;
            }
            rfIdx++;
        }
        else  //rfIdx == vPos
        {
            //search buffered text
            segment.clear();
            if (tBuff.length() > 0) {
                segment.push_back(tBuff);
                multiedsm->searchNextSegment(segment);
                tBuff = "";
            }

            //then search current variant segment and skip required number of characters
            multiedsm->searchNextSegment(vSeg);
            unsigned int j;
            for (j = 1; j < vSkip; j++) {
                getNextCharFAS(rf);
            }
            rfIdx += vSkip;

            //get the next variant position, skipping repeats or nested variants
            do {
                if (++vit == variantItems.end()) {
                    break;
                }
                vPos = vit->pos;
                vSkip = vit->skip;
                vSeg = vit->seg;
            }
            while (vPos < rfIdx);
        }
    }
    if (tBuff.length() > 0)
    {
        segment.clear();
        segment.push_back(tBuff);
        multiedsm->searchNextSegment(segment);
        tBuff = "";
    }

    rf.close();

    return true;
}

/**
 * Search for patterns using the EDS (Elastic Degenerate Sequence) format file
 *
 * @param multiedsm The multiedsm object (reference)
 * @param seqFile The EDS sequence file name
 * @return 0=Success, 1=Error
 */
int searchEDS(MultiEDSM * multiedsm, string & seqFile)
{
    //open sequence file
    ifstream eds(seqFile.c_str(), ios::in);
    if (!eds.good()) {
        cerr << "Error. Unable to open sequence file!" << endl;
        return false;
    }

    //initialize variables required for searching
    Segment tempSeg;
    string x = "";
    x.reserve(BUFFERSIZE);
    char c = 0;
    bool inDegSeg = false;
    unsigned int i = 0, j = 0;

    //go through the sequence file
    while ((c = getNextCharEDS(eds)) != '\0')
    {
        if (i == 0 && c == '{')
        {
            inDegSeg = true;
        }
        else if (c == ',')
        {
            tempSeg.push_back(x);
            x = "";
            j = 0;
        }
        else if (c == '}' || (c == '{' && i > 0))
        {
            if (x.length() > 0) {
                tempSeg.push_back(x);
                x = "";
                j = 0;
                multiedsm->searchNextSegment(tempSeg);
                tempSeg.clear();
            }
            inDegSeg = (c == '{');
        }
        else if (c != '{')
        {
            x += c;
            j++;
            if (!inDegSeg && j == BUFFERSIZE)
            {
                tempSeg.push_back(x);
                x = "";
                j = 0;
                multiedsm->searchNextSegment(tempSeg);
                tempSeg.clear();
            }
        }
        i++;
    }
    if (x != "") {
        tempSeg.push_back(x);
        multiedsm->searchNextSegment(tempSeg);
        x = "";
        tempSeg.clear();
    }

    eds.close();

    return true;
}

/**
 * Main.
 */
int main(int argc, char * argv[])
{
    //parse command-line options
    int c, x = 0, optind = 1;
    EDSDEGLENTYPE countingType;
    unsigned int stfL = 1;
    string seqF, varF, patF, memL, cntL;
    static struct option long_options[] = {
        {"sequence-file", required_argument, 0, 's'},
        {"variants-file", required_argument, 0, 'v'},
        {"patterns-file", required_argument, 0, 'p'},
        {"mem-limit",     required_argument, 0, 'm'},
        {"stf-limit",     required_argument, 0, 't'},
        {"counting-type", required_argument, 0, 'c'},
        {"help",          no_argument,       0, 'h'},
        {0,               0,                 0,  0 }
    };

    string help = "Multiple Elastic Degenerate String Matching (Multi-EDSM) ---\n\n\
Example EDS Type Search: ./multiedsm -s seq.eds -p patterns.txt -m 4g\n\
Example FASTA+VCF Type Search: ./multiedsm -s reference.fasta -v variants.vcf -p patterns.txt -m 4g\n\n\
Standard (Required) Arguments:\n\
  -s\t--sequence-file\t<str>\tThe EDS or reference FASTA file. Use the correct file for the correct search type.\n\
  -v\t--variants-file\t<str>\tThe VCF variants-file. Support for .vcf.gz. Use only with FASTA+VCF search type.\n\
  -p\t--patterns-file\t<str>\tThe patterns file. Each pattern must be on a different line.\n\
  -m\t--mem-limit\t<str>\tThe maximum amount of memory to use. Use 'g' or 'm' modifiers for GB or MB, e.g. 3.5g\n\n\
Optional Arguments:\n\
  -t\t--stf-limit\t<uint>\tManually set the Suffix Tree-based data structure memory usage-limiting factor -  \n\
    \t           \t      \tset it to 1 to use as much memory as possible, O to use O(M) memory, or anything up\n\
    \t           \t      \tto O([M/w]) for balanced memory usage and speed. (default=1)\n\
  -c\t--counting-type\t<str>\tCount degenerate segments as one position (FIXEDLENGTH),\n\
    \t               \t     \tthe length of the first string (FIRSTLENGTH), or\n\
    \t               \t     \tup to the length of the first string (UPTOFIRSTLENGTH=default)?\n\n\
Miscellaneous:\n\
  -h\t--help\t\t<void>\tThis help message.\n";

    while ((c = getopt_long(argc, argv, "s:v:p:m:t:c:h", long_options, &optind)) != -1)
    {
        switch (c)
        {
            case 's':
                seqF = optarg;
                x++;
                break;
            case 'v':
                varF = optarg;
                x++;
                break;
            case 'p':
                patF = optarg;
                x++;
                break;
            case 'm':
                memL = optarg;
                x++;
                break;
            case 't':
                stfL = (unsigned int) atoi(optarg);
                break;
            case 'c':
                cntL = optarg;
                break;
            case 'h':
                cout << help << endl;
                return EXIT_SUCCESS;
            default:
                cerr << "Error: unrecognised argument." << endl << help << endl;
                return EXIT_FAILURE;
        }
    }

    if (x == 3)
    {
        if (seqF == "" || patF == "" || memL == "") {
            cerr << "Error: Invalid arguments!" << endl << help << endl;
            return EXIT_FAILURE;
        }
    }
    else if (x == 4)
    {
        if (seqF == "" || varF == "" || patF == "" || memL == "") {
            cerr << "Error: Invalid arguments!" << endl << help << endl;
            return EXIT_FAILURE;
        }
    }
    else
    {
        cerr << "Error: Invalid number of arguments!" << endl << help << endl;
        return EXIT_FAILURE;
    }

    if (cntL == "FIXEDLENGTH") {
        countingType = EDSDEGLENTYPE::FIXEDLENGTH;
    } else if (cntL == "FIRSTLENGTH") {
        countingType = EDSDEGLENTYPE::FIRSTLENGTH;
    } else if (cntL == "UPTOFIRSTLENGTH" || cntL == "") {
        countingType = EDSDEGLENTYPE::UPTOFIRSTLENGTH;
    } else {
        cerr << "Error: Invalid counting type given!" << endl << help << endl;
        return EXIT_FAILURE;
    }

    unsigned long long int memLimit = 1024 * 1024 * (unsigned long long int)getMemLimitMB(memL); //memLimit in bytes
    if (memLimit == 0) {
        cerr << "Error: Invalid memory value given!" << endl << help << endl;
        return EXIT_FAILURE;
    }

    //open patterns file and determine total size of patterns in order to calculate memory requirements
    unsigned long long int M = 0, k = 0, maxP = 0;
    ifstream pf(patF.c_str(), ios::in);
    if (!pf.good()) {
        cerr << "Error: Failed to open pattern file!" << endl;
        return EXIT_FAILURE;
    }
    vector<string> patterns;
    string pattern;
    while(getline(pf, pattern))
    {
        if (pattern.length() > maxP) {
            maxP = pattern.length();
        }
        M += pattern.length();
        patterns.push_back(pattern);
        k++;
        if (M + k >= memLimit) {
            pf.close();
            cerr << "Error: Set of patterns too large -- memory limit exceeded! Try splitting patterns up into batches or increasing the memory limit." << endl;
            return EXIT_FAILURE;
        }
    }
    pf.close();

    //calculate base memory requirements making sure there is enough memory for at least two levels of the OccVector
    unsigned long long int i, stpnodes, suffixtree, stp2pos, pos2pat, shiftand, twentyonebitvectors, total;
    stpnodes = M + k;                                                                    //pattern string for STp
    suffixtree = stpnodes * 24;                                                          //STp - num nodes < 2n, suffix array 6n, so 2n * 6n = ~12n and same for suffix tree and associated datastructures so ~24n altogether
    stp2pos = sizeof(unsigned int) * stpnodes;                                           //stp2pos
    pos2pat = sizeof(unsigned int) * M;                                                  //pattern positions vector
    shiftand = (unsigned long long int) ceil((double)M / (double)BITSINWORD) * (3 + SIGMA) * WORDSIZE;   //shiftand, sigma & Sv, Ev and D
    twentyonebitvectors = (unsigned long long int) ceil((double)M / (double)BITSINWORD) * WORDSIZE * 21; //level1 = 4 BVs, level2 = 17 BVs, level1 + level2 = 21 BVs
    total = M + stpnodes + suffixtree + stp2pos + pos2pat + shiftand + twentyonebitvectors + BUFFERSIZE;
    if (total >= memLimit) {
        cerr << "Error: Insufficient memory to continue! Try increasing the memory limit." << endl;
        return EXIT_FAILURE;
    }

    //calculate left-over memory for use in raw bitvector storage
    total = total - twentyonebitvectors;
    unsigned long long int remainingMemory = memLimit - total;
    unsigned long long int maxNoBitVectorsStorable = (unsigned long long int) ceil((double)remainingMemory / (ceil((double)M / (double)BITSINWORD) * WORDSIZE));

    //start MultiEDSM search
    cout << "Multi-EDSM starting..." << endl << endl;
    MultiEDSM * multiedsm = new MultiEDSM(ALPHABET, patterns, maxNoBitVectorsStorable, stfL, countingType);
    patterns.clear();

    cout << "Searching..." << endl;

    bool success;
    if (x == 4) {
        success = searchVCF(multiedsm, seqF, varF);
    } else {
        success = searchEDS(multiedsm, seqF);
    }
    if (!success) {
        cerr << "Error: An error occured." << endl;
        return EXIT_FAILURE;
    }

    cout << "No. determinate bases (f): " << multiedsm->getf() << endl;
    cout << "No. degenerate bases (F): " << multiedsm->getF() << endl;
    cout << "No. determinate segments (d): " << multiedsm->getd() << endl;
    cout << "No. degenerate segments (D): " << multiedsm->getD() << endl;
    cout << "No. strings processed shorter than pattern (N'): " << multiedsm->getNp() << endl;
    cout << "Multi-EDSM total processing time: " << multiedsm->getDuration() << "s." << endl << endl;

    if (multiedsm->getMatches().size() == 0)
    {
        cout << "No matches found." << endl;
    }
    else
    {
        cout << multiedsm->getMatches().size() << " matches found!" << endl << endl;
        cout << "Position,PatternId" << endl << "------------------" << endl;
        for (const pair<unsigned int, unsigned int> & match : multiedsm->getMatches()) {
            cout << match.first << "," << match.second << endl;
        }
    }

    delete multiedsm;

    return EXIT_SUCCESS;
}
