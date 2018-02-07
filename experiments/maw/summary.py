#
#    Copyright (C) 2017 Ahmad Retha and Solon P. Pissis.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import os

MAWFile = './patterns/filtered_hg37.maws'
originalMAWs = []
with open(MAWFile, 'r') as f:
    for motif in f:
        originalMAWs.append(motif.strip())

print 'No. MAWs originally extracted: ' + str(len(originalMAWs))

chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y']
falseMAWsFile = './results/falseMAWs_%s.txt'
for chromosome in chromosomes:
    chrFile = falseMAWsFile % chromosome
    if os.path.exists(chrFile):
        with open(chrFile, 'r') as f:
            for index in f:
                originalMAWs[int(index)] = '-'
    else:
        print 'Warning: Missing file ' + chrFile

confirmedMAWs = set(originalMAWs)
confirmedMAWs.remove('-')

confirmedMAWsFile = './results/confirmedMAWs.maws'
with open(confirmedMAWsFile, 'w') as f:
    j = 0
    for maw in confirmedMAWs:
        if j != 0:
            f.write('\n')
        f.write(maw)
        j += 1

print 'No. MAWs confirmed: ' + str(j) + ' (%0.2f%%)' % (float(j)/float(len(originalMAWs)) * 100.0)
print 'Confirmed MAWs written to: ' + confirmedMAWsFile
