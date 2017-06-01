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
 * Convert raw string like 3.5g (3.5 gigabytes) to 3584 (megabytes)
 */
unsigned int getMemLimitMB(string rawString)
{
    int l = rawString.length();
    if (l == 0) {
        return UINT_MAX;
    }
    char unit = rawString[l - 1];
    if (unit == 'm' || unit == 'M') {
        return (unsigned int) ceil(atof(rawString.substr(0, l - 1).c_str()));
    } else if (unit == 'g' || unit == 'G') {
        return (unsigned int) ceil(atof(rawString.substr(0, l - 1).c_str()) * 1024);
    } else {
        cerr << "Invalid memory argument used: \"" << rawString << "\"! Continuing with unlimited memory." << endl;
        return UINT_MAX;
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
  -p\t--patterns-file\t<str>\tThe patterns file. Each pattern must be on a different line.\n\n\
Optional Arguments:\n\
  -m\t--mem-limit\t<float>\tThe maximum amount of memory to use. Use 'g' or 'm' modifiers, e.g. 3.5g\n\n\
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
                break;
            case 'h':
                cout << help << endl;
                return 0;
            default:
                cerr << "Error: unrecognised argument." << endl << help << endl;
                return 1;
        }
    }

    if (x == 2)
    {
        if (seqF == "" || patF == "") {
            cerr << "Error: Invalid arguments!" << endl << help << endl;
            return 1;
        }
    }
    else if (x == 3)
    {
        if (seqF == "" || varF == "" || patF == "") {
            cerr << "Error: Invalid arguments!" << endl << help << endl;
            return 1;
        }
    }
    else
    {
        cerr << "Error: Invalid number of arguments!" << endl << help << endl;
        return 1;
    }

    //open patterns file and search patterns in batches
    bool success;
    unsigned int M = 0, k = 0, batchNo = 1, batchCount = 0, memLimit = getMemLimitMB(memL);
    unsigned int f = 0, F = 0, d = 0, D = 0, Np = 0;
    double duration = 0.0;
    unsigned int stp, stp2pos, pos2pat, ovmem, shiftand, total;
    vector<pair<unsigned int, unsigned int>> foundList;
    vector<string> patterns;
    string pattern;
    ifstream pf(patF.c_str(), ios::in);
    if (!pf.good()) {
        cerr << "Error: Failed to open pattern file!" << endl;
        return 1;
    }
    cout << "Multi-EDSM started..." << endl << endl;
    while(getline(pf, pattern))
    {
        M += pattern.length();
        k++;
        stp = M + k;
        stp2pos = sizeof(unsigned int) * stp;
        pos2pat = sizeof(unsigned int) * M;
        ovmem = (unsigned int) ceil(stp / BITSINWORD) * stp * WORDSIZE;
        shiftand = (unsigned int) ceil(M / BITSINWORD) * (3 + SIGMA) * WORDSIZE;
        patterns.push_back(pattern);
        total = (unsigned int) ((M + stp + stp2pos + pos2pat + ovmem + shiftand + BUFFERSIZE) / (1024 * 1024));
        if (total >= memLimit)
        {
            //perform the search
            MultiEDSM multiedsm(ALPHABET, patterns);
            if (x == 3) {
                success = searchVCF(&multiedsm, seqF, varF);
            } else {
                success = searchEDS(&multiedsm, seqF);
            }
            if (!success) {
                return 1;
            }
            cout << "Searched batch " << batchNo << " containing " << k << " patterns." << endl;
            //grab the results and update statistics
            for (const pair<unsigned int, unsigned int> & match : multiedsm.getMatches()) {
                foundList.push_back(pair<unsigned int, unsigned int>(match.first, match.second + batchCount));
            }
            f += multiedsm.getf();
            F += multiedsm.getF();
            d += multiedsm.getd();
            D += multiedsm.getD();
            Np += multiedsm.getNp();
            duration += multiedsm.getDuration();
            //reset batch searching variables
            M = 0;
            k = 0;
            batchNo++;
            batchCount += patterns.size();
            patterns.clear();
        }
    }
    pf.close();

    //search the remaining patterns
    MultiEDSM multiedsm(ALPHABET, patterns);
    if (patterns.size() > 0)
    {
        //perform the search
        if (x == 3) {
            success = searchVCF(&multiedsm, seqF, varF);
        } else {
            success = searchEDS(&multiedsm, seqF);
        }

        if (!success) {
            return 1;
        }
        if (batchNo > 1) {
            cout << "Searched batch " << batchNo << " containing " << k << " patterns." << endl;
        }
        //grab the results and update statistics
        for (const pair<unsigned int, unsigned int> & match : multiedsm.getMatches()) {
            foundList.push_back(pair<unsigned int, unsigned int>(match.first, match.second + batchCount));
        }
        f += multiedsm.getf();
        F += multiedsm.getF();
        d += multiedsm.getd();
        D += multiedsm.getD();
        Np += multiedsm.getNp();
        duration += multiedsm.getDuration();
        //reset batch searching variables
        patterns.clear();
        M = 0;
        k = 0;
    }

    //output results
    cout << endl;
    if (batchNo > 1)
    {
        cout << "Patterns searched in " << batchNo << " batches." << endl;
        cout << "No. determinate bases (f): (" << f << " processed) " << multiedsm.getf() << endl;
        cout << "No. degenerate bases (F): (" << F << " processed) " << multiedsm.getF() << endl;
        cout << "No. determinate segments (d): (" << d << " processed) " << multiedsm.getd() << endl;
        cout << "No. degenerate segments (D): (" << D << " processed) " << multiedsm.getD() << endl;
    }
    else
    {
        cout << "No. determinate bases (f): " << f << endl;
        cout << "No. degenerate bases (F): " << F << endl;
        cout << "No. determinate segments (d): " << d << endl;
        cout << "No. degenerate segments (D): " << D << endl;
    }
    cout << "No. strings processed shorter than pattern (N'): " << Np << endl;
    cout << "EDSM-BV processing time: " << duration << "s." << endl << endl;

    if (foundList.size() == 0)
    {
        cout << "No matches found." << endl;
    }
    else
    {
        cout << foundList.size() << " matches found!" << endl << endl;
        cout << "Position,PatternId" << endl << "------------------" << endl;
        for (const pair<unsigned int, unsigned int> & match : foundList) {
            cout << match.first << "," << match.second << endl;
        }
    }

    return 0;
}
