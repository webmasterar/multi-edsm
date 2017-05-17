
import os
import sys
import random
import argparse
import subprocess

#
# Parse command line arguments and check input
#
parser = argparse.ArgumentParser(description='Test Multi-EDSM using randomly-generated synthetic data. Provide \
a pattern size like 16 and the number of patterns you want to randomly generate and it will run Multi-EDSM on \
the provided sequence file.')
parser.add_argument('pattern_size', type=int, help='the length of each randomly generated pattern, e.g. 16')
parser.add_argument('total_number', type=int, help='How many patterns to generate, e.g. 100')
parser.add_argument('sequence_file', type=argparse.FileType('r'), help='The path to the sequence file to search through')
parser.add_argument('--recreate-files', action='store_true', help='Delete and recreate the pattern file')
args = parser.parse_args()

if args.total_number < 1:
	print 'Error: You specified an invalid number of patterns to generate!'
	sys.exit(1)

sequenceFileName = os.path.basename(args.sequence_file.name)
sequenceFileName = ".".join(sequenceFileName.split('.')[0:-1])
patternsFile = './randomPatterns/' + str(args.pattern_size) + '_' + str(args.total_number) + '.txt'
resultsFile = './results/' + str(args.pattern_size) + '_' + str(args.total_number) + '_' + sequenceFileName + '.txt'

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
	count = 0
	with open(patternsFile, 'w') as f:
		while count < args.total_number:
			if count != 0:
				f.write("\n")
			for j in range(args.pattern_size):
				f.write(alphabet[random.randint(0, sigma - 1)])
			count += 1

#
# Search for the patterns with Multi-EDSM
#
cmd = './../../multiedsm ' + args.sequence_file.name + ' ' + patternsFile
output = subprocess.Popen(['time ' + cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
with open(resultsFile, 'w') as f:
	f.write(output)

print output
