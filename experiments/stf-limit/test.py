
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
