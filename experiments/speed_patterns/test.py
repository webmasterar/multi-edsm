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
import random
import argparse
import subprocess

#
# Parse command line arguments and check input
#
parser = argparse.ArgumentParser(description='Test Multi-EDSM using randomly-generated synthetic data. Provide \
a pattern size like 16 and a total length like 64 and it will create 4 random patterns with combined length 64 \
and run Multi-EDSM on the provided sequence file.')
parser.add_argument('pattern_size', type=int, help='the length of each randomly generated pattern, e.g. 16')
parser.add_argument('num_patterns', type=int, help='How many patterns to make, e.g. 100')
parser.add_argument('sequence_file', type=argparse.FileType('r'), help='The path to the sequence file to search through')
parser.add_argument('--recreate-files', action='store_true', help='Delete and recreate the pattern file')
args = parser.parse_args()

sequenceFileName = os.path.basename(args.sequence_file.name)
sequenceFileName = ".".join(sequenceFileName.split('.')[0:-1])
patternsFile = './randomPatterns/' + str(args.pattern_size) + '_' + str(args.num_patterns) + '.txt'
resultsFile = './results/' + str(args.pattern_size) + '_' + str(args.num_patterns) + '_' + sequenceFileName + '.txt'

if args.recreate_files:
	if os.path.isfile(patternsFile):
		os.remove(patternsFile)
	if os.path.isfile(resultsFile):
		os.remove(resultsFile)

#
# Create the patterns file
#
alphabet = 'ACGT'
sigma = len(alphabet)
if not os.path.exists(patternsFile):
	i = 0
	with open(patternsFile, 'w') as f:
		while i < args.num_patterns:
			if i != 0:
				f.write("\n")
			for j in range(args.pattern_size):
				f.write(alphabet[random.randint(0, sigma - 1)])
			i += 1

#
# Search for the patterns with Multi-EDSM
#
cmd = './../../multiedsm -m 4g -s ' + args.sequence_file.name + ' -p ' + patternsFile
output = subprocess.Popen(['time ' + cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
with open(resultsFile, 'w') as f:
	f.write(output)

print output
