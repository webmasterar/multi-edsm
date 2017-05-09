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

#include <vector>
#include <string>
#include "UnrestrictedMultiShiftAnd.hpp"

/**
 * MyUMSA is a subclass of the UnrestrictedMultiShiftAnd class and it adds an
 * extra function: getEndingStates()
 */
class MyUMSA : public UnrestrictedMultiShiftAnd
{
public:
    /**
     * @constructor Calls superclass -- see UnrestrictedMultiShiftAnd doc for more info
     * @param alphabet
     */
    MyUMSA(const std::string & alphabet) : UnrestrictedMultiShiftAnd(alphabet) {};

    /**
     * This function simply returns a WordVector containing the ending states of
     * the Shift-And patterns, where the 1s are indicators of matching positions
     * for when a pattern successfully matches in a searched text
     *
     * @return A WordVector marking the ending positions/states of the pattern set
     */
    const std::vector<WORD> getEndingStates() const
    {
        return this->Ev;
    }
};

#endif
