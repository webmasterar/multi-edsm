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

char BUFF[BUFFERSIZE];
int BUFFLIMIT = 0;
int POS = 0;

/*
* Buffered file reading char by char
*
* @param f opened file handle
*/
char getNextChar(ifstream & f)
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
    return c;
}

/**
 * Convert raw string like 3.5g (3.5 gigabytes) to 3584 (megabytes) -- return 0 indicated error
 */
unsigned int getMemLimitMB(string rawString)
{
    int l = rawString.length();
    if (l == 0) {
        return EXIT_SUCCESS;
    }
    char unit = rawString[l - 1];
    if (unit == 'm' || unit == 'M') {
        return (unsigned int) ceil(atof(rawString.substr(0, l - 1).c_str()));
    } else if (unit == 'g' || unit == 'G') {
        return (unsigned int) ceil(atof(rawString.substr(0, l - 1).c_str()) * 1024);
    } else {
        return EXIT_SUCCESS;
    }
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
    ifstream rf(seqFile.c_str(), ios::in);
    if (!rf.good()) {
        cerr << "Error: Failed to open reference file!" << endl;
        return false;
    }

    VariantCallFile vf;
    vf.open(varFile);
    if (!vf.is_open()) {
        cerr << "Error: Failed to open variants file!" << endl;
        return false;
    }

    //initialize fasta file reading and segment creation helper variables
    string tBuff = "";
    tBuff.reserve(BUFFERSIZE);
    char c;
    unsigned int rfIdx = 1, vfIdx = 0, i = 0;
    Segment segment;

    //skip first line of fasta file
    getline(rf, tBuff);
    tBuff = "";

    //create variables for reading through vcf records and looking for duplicates
    Variant var(vf), vBuffer(vf), var2(vf);
    bool hasMoreVariants = true;
    Segment vAlleles;

    //read first variant and possibly successive duplicates for the same position, removing duplicate alleles
    hasMoreVariants = vf.getNextVariant(var);
    if (hasMoreVariants)
    {
        vfIdx = (unsigned int) var.position;
        for (const auto & a : var.alleles) {
            if (a[0] != '<') {
                vAlleles.push_back(a);
            }
        }
        while (true) {
            hasMoreVariants = vf.getNextVariant(var2);
            if (hasMoreVariants)
            {
                if (var2.position == var.position) {
                    for (const auto & a : var2.alt) {
                        if (a[0] != '<') {
                            vAlleles.push_back(a);
                        }
                    }
                } else {
                    vBuffer = var2;
                    break;
                }
            }
            else
            {
                break;
            }
        }
    }

    //go through the reference sequence
    while ((c = getNextChar(rf)) != '\0')
    {
        if (!(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')) {
            continue;
        }

        if (rfIdx != vfIdx)
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
        }
        else
        {
            segment.clear();
            if (tBuff.length() > 0) {
                segment.push_back(tBuff);
                multiedsm->searchNextSegment(segment);
                tBuff = "";
            }

            //then search current variant
            if (vAlleles.size() > 0) {
                multiedsm->searchNextSegment(vAlleles);
                vAlleles.clear();
            }

            //fetch the next variant to be searched for when its position comes up
            if (vBuffer.alleles.size() > 0)
            {
                vfIdx = (unsigned int) vBuffer.position;
                for (const auto & a : vBuffer.alleles) {
                    if (a[0] != '<') {
                        vAlleles.push_back(a);
                    }
                }

                vBuffer.alleles.clear();

                while (true) {
                    hasMoreVariants = vf.getNextVariant(var2);
                    if (hasMoreVariants)
                    {
                        if (vfIdx == (unsigned int) var2.position) {
                            for (const auto & a : var2.alt) {
                                if (a[0] != '<') {
                                    vAlleles.push_back(a);
                                }
                            }
                        } else {
                            vBuffer = var2;
                            break;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        rfIdx++;
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
    while ((c = getNextChar(eds)) != '\0')
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
            switch (c) {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'N':
                case EPSILON[0]:
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
                    break;
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
    string seqF, varF, patF, memL;
    static struct option long_options[] = {
        {"sequence-file", required_argument, 0, 's'},
        {"variants-file", required_argument, 0, 'v'},
        {"patterns-file", required_argument, 0, 'p'},
        {"mem-limit",     required_argument, 0, 'm'},
        {"help",          no_argument,       0, 'h'},
        {0,               0,                 0, 0}
    };

    string help = "Multiple Elastic Degenerate String Matching (Multi-EDSM) ---\n\n\
Example EDS Type Search: ./multiedsm --sequence-file seq.eds --patterns-file patterns.txt\n\
Example FASTA+VCF Type Search: ./multiedsm --sequence-file reference.fasta --variants-file variants.vcf --patterns-file patterns.txt\n\n\
Standard (Required) Arguments:\n\
  -s\t--sequence-file\t<str>\tThe EDS or reference FASTA file. Use the correct file for the correct search type.\n\
  -v\t--variants-file\t<str>\tThe VCF variants-file. Support for .vcf.gz. Use only with FASTA+VCF search type.\n\
  -p\t--patterns-file\t<str>\tThe patterns file. Each pattern must be on a different line.\n\
  -m\t--mem-limit\t<str>\tThe maximum amount of memory to use. Use 'g' or 'm' modifiers for GB or MB, e.g. 3.5g\n\n\
Miscellaneous:\n\
  -h\t--help\t\t<void>\tThis help message.\n";

    while ((c = getopt_long(argc, argv, "s:v:p:m:h", long_options, &optind)) != -1)
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
    // cout << "stpnodes:" << stpnodes << " suffixtree:" << suffixtree << " stp2pos:" << stp2pos << " pos2pat:" << pos2pat << " shiftand:" << shiftand << " twentyonebitvectors:" << twentyonebitvectors << " total:" << total << endl;
    if (total >= memLimit) {
        cerr << "Error: Insufficient memory to continue! Try increasing the memory limit." << endl;
        return EXIT_FAILURE;
    }

    //calculate left-over memory for use in raw bitvector storage
    total = total - twentyonebitvectors;
    unsigned long long int remainingMemory = memLimit - total;
    unsigned long long int maxNoBitVectorsStorable = (unsigned long long int) ceil((double)remainingMemory / (ceil((double)M / (double)BITSINWORD) * WORDSIZE));
    // cout << "Memlimit:" << memLimit << " and remainingMemory:" << remainingMemory << " maxNoBitVectorsStorable:" << maxNoBitVectorsStorable << endl;

    //start MultiEDSM search
    cout << "Multi-EDSM starting..." << endl << endl;
    MultiEDSM * multiedsm = new MultiEDSM(ALPHABET, patterns, maxNoBitVectorsStorable);
    patterns.clear();
    bool success;
    if (x == 4) {
        success = searchVCF(multiedsm, seqF, varF);
    } else {
        success = searchEDS(multiedsm, seqF);
    }
    if (!success) {
        cerr << "Error: An unspecified error occured." << endl;
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
