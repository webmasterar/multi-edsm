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

int main(int argc, char * argv[])
{
    string help = "There are two ways to run Multiple Elastic Degenerate String Matching (Multi-EDSM) ---\n\
    \tUsage: ./multiedsm seq.txt patterns.txt\n\
    \tUsage: ./multiedsm reference.fasta variants.vcf patterns.txt";

    if (argc == 1 || (argc == 2 && (strcmp("--help", argv[1]) == 0 || strcmp("-h", argv[1]) == 0))) {
        cout << help << endl;
        return 0;
    }

    if (!(argc == 4 || argc == 3)) {
        cerr << "Invalid number of arguments!" << endl;
        cout << help << endl;
        return 1;
    }

    //pattern file
    string pattFile;
    if (argc == 3) {
        pattFile = argv[2];
    } else {
        pattFile = argv[3];
    }
    ifstream pf(pattFile.c_str(), ios::in);
    if (!pf.good()) {
        cerr << "Error: Failed to open pattern file!" << endl;
        return 1;
    }

    //patterns
    vector<string> patterns;
    string pattern;
    while(getline(pf, pattern)) {
        patterns.push_back(pattern);
    }
    pf.close();

    //MultiEDSM
    MultiEDSM multiedsm("ACGT", patterns);
    cout << "Multi-EDSM searching..." << endl << endl;

    //Searching Reference+VCF file
    if (argc == 4)
    {
        string refName = argv[1];
        ifstream rf(refName.c_str(), ios::in);
        if (!rf.good()) {
            cerr << "Error: Failed to open reference file!" << endl;
            return 1;
        }

        string vcfName = argv[2];
        VariantCallFile vf;
        vf.open(vcfName);
        if (!vf.is_open()) {
            cerr << "Error: Failed to open variants file!" << endl;
            return 1;
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
                    multiedsm.searchNextSegment(segment);
                    tBuff = "";
                    i = 0;
                }
            }
            else
            {
                segment.clear();
                if (tBuff.length() > 0) {
                    segment.push_back(tBuff);
                    multiedsm.searchNextSegment(segment);
                    tBuff = "";
                }

                //then search current variant
                if (vAlleles.size() > 0) {
                    multiedsm.searchNextSegment(vAlleles);
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
            multiedsm.searchNextSegment(segment);
            tBuff = "";
        }

        rf.close();
    }
    //searching custom EDS format file
    else
    {
        //open sequence file
        ifstream eds(argv[1], ios::in);
        if (!eds.good()) {
            cerr << "Error. Unable to open sequence file!" << endl;
            return 1;
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
                    multiedsm.searchNextSegment(tempSeg);
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
                            multiedsm.searchNextSegment(tempSeg);
                            tempSeg.clear();
                        }
                        break;
                }
            }
            i++;
        }
        if (x != "") {
            tempSeg.push_back(x);
            multiedsm.searchNextSegment(tempSeg);
            x = "";
            tempSeg.clear();
        }

        eds.close();
    }

    //output results
    cout << endl;
    cout << "No. determinate bases (f): " << multiedsm.getf() << endl;
    cout << "No. degenerate bases (F): " << multiedsm.getF() << endl;
    cout << "No. determinate segments (d): " << multiedsm.getd() << endl;
    cout << "No. degenerate segments (D): " << multiedsm.getD() << endl;
    cout << "No. strings processed shorter than pattern (N'): " << multiedsm.getNp() << endl;
    cout << "EDSM-BV processing time: " << multiedsm.getDuration() << "s." << endl << endl;

    if (multiedsm.getMatches().size() == 0)
    {
        cout << "No matches found." << endl;
    }
    else
    {
        cout << multiedsm.getMatches().size() << " matches found!" << endl << endl;
        cout << "Position,PatternId" << endl << "------------------" << endl;
        for (const pair<unsigned int, unsigned int> & match : multiedsm.getMatches()) {
            cout << match.first << "," << match.second << endl;
        }
    }

    return 0;
}
