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
import sys
import math
import argparse
import subprocess

#
# Parse command line arguments and check input
#
parser = argparse.ArgumentParser(description='Measure performance of Multi-EDSM with changing SuffixTree Factor limit up to [M/w]')
parser.add_argument('patternsFile', type=argparse.FileType('r'), help='The patterns file')
parser.add_argument('edsFile', type=argparse.FileType('r'), help='The EDS file')
parser.add_argument('steps', type=int, default=5, help='The number of steps to test up to [M/w] (default = 5)')
parser.add_argument('--memory', type=int, default=4096, help='How much RAM to use in MB (default = 4096)')
args = parser.parse_args()

w = 64
M = 0
for pattern in args.patternsFile:
    M += len(pattern.strip())
step = int(math.ceil( math.ceil(float(M)/float(w)) / float(args.steps) ))

print 'Running experiment in %d steps for M=%d' % (args.steps, M)
print ''

#
# Search for the patterns with Multi-EDSM and output results
#
for i in range(1, args.steps+1):
    print 'Step %d: stf-limit=%d' % (i, i * step)
    cmd = './../../multiedsm -m ' + str(args.memory) + 'm -s ' + args.edsFile.name + ' -p ' + args.patternsFile.name + ' -t ' + str(i * step)
    print cmd
    output = subprocess.Popen(['time ' + cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    print output

print 'Done!'
