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
#include "MultiEDSM.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    vector<string> patterns;
    patterns.push_back("ACACA");
    patterns.push_back("CACCA");
    MultiEDSM multiedsm("ACGT", patterns);
    vector<string> segment;
    segment.push_back("CAACACA");
    segment.push_back("E");
    segment.push_back("CACC");
    multiedsm.searchNextSegment(segment);
    segment.clear();
    segment.push_back("AA");
    segment.push_back("AC");
    multiedsm.searchNextSegment(segment);
    segment.clear();
    segment.push_back("CA");
    multiedsm.searchNextSegment(segment);
    return 0;
}
