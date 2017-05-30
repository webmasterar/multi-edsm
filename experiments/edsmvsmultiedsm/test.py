
import os
import sys
import time
import random
import argparse
import subprocess

#
# Parse command line arguments and check input
#
parser = argparse.ArgumentParser(description='Test Multi-EDSM\'s speed against \
EDSM-BV by searching many randomly-generated patterns. Tell test.py which \
algorithm to run, multi=Multi-EDSM or single=EDSM-BV, the pattern size to search \
for e.g. 40 and the number of patterns you wish to test both programs with.')
parser.add_argument('algorithm', choices=['multi','single'])
parser.add_argument('pattern_size', type=int, help='the length of each randomly generated pattern, e.g. 40')
parser.add_argument('total_number', type=int, help='How many patterns to search for, e.g. 500')
parser.add_argument('sequence_file', type=argparse.FileType('r'), help='The path to the sequence file to search through')
parser.add_argument('--recreate-files', action='store_true', help='Delete and recreate the pattern file')
args = parser.parse_args()

if args.pattern_size < 2 or args.pattern_size > 63:
	print 'Error: Invalid pattern size'
	sys.exit(1)
if args.total_number < 1:
	print 'Error: Invalid number of patterns'
	sys.exit(1)

sequenceFileName = os.path.basename(args.sequence_file.name)
sequenceFileName = ".".join(sequenceFileName.split('.')[0:-1])
patternsFile = './randomPatterns/' + str(args.pattern_size) + '_' + str(args.total_number) + '.txt'
resultsFile = './results/' + args.algorithm + '_' + str(args.pattern_size) + '_' + str(args.total_number) + '_' + sequenceFileName + '.txt'

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
	total = 0
	with open(patternsFile, 'w') as f:
		while total < args.total_number:
			if total != 0:
				f.write("\n")
			for j in range(args.pattern_size):
				f.write(alphabet[random.randint(0, sigma - 1)])
			total += 1

#
# Search for the patterns
#
if args.algorithm == 'multi':
	#
	# Search for the patterns with Multi-EDSM
	#
	cmd = './../../multiedsm ' + args.sequence_file.name + ' ' + patternsFile
	output = subprocess.Popen(['time ' + cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
	with open(resultsFile, 'w') as f:
		f.write(output)
	print output

else:
	#
	# Search for the patterns with EDSM-BV
	#
	with open(resultsFile, 'w') as r:
		with open(patternsFile, 'r') as f:
			totalTime = 0.0
			for pattern in f:
				cmd = './../../../edsm/edsm ' + args.sequence_file.name + ' ' + pattern
				startTime = time.time()
				try:
					output = subprocess.Popen(['time ' + cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
				finally:
					endTime = time.time()
				totalTime += endTime - startTime
				output += '\nEDSM-BV running time: ' + str(endTime - startTime) + '\n\n'
				r.write(output)
				print output
			tot = '\n\nTotal Running Time of EDSM-BV: ' + str(totalTime)
			r.write(tot)
			print tot
