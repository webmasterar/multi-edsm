
import os
import sys
import random
import argparse
import subprocess

#
# Takes a raw bit of EDS sequence such as "A,GG}CACA{G,GG}A" and returns a non-
# degenerate pattern that can be used in a search, e.g. ACACAGGA.
#
def processRawPattern(rawPattern):
	p = ''
	l = len(rawPattern)
	i = 0
	while i < l:
		if rawPattern[i] == '{':
			i += 1
			degSeg = ''
			while i < l and rawPattern[i] != '}':
				degSeg += rawPattern[i]
				i += 1
			degSegArr = degSeg.split(',')
			p += degSegArr[random.randint(0, len(degSegArr) - 1)]
		elif rawPattern[i] == ',':
			while i < l and rawPattern[i] != '}':
				i += 1
		elif rawPattern[i] != '}':
			p += rawPattern[i]
		i += 1
	return p

#
# Parse command line arguments and check input
#
parser = argparse.ArgumentParser(description='Test Multi-EDSM on real data. Provide \
the input EDS format file and how many patterns you wish to test it with -- \
patterns are automatically extracted from the EDS file. Optionally, specify the \
minimum and maximum pattern length with the --min-size and --max-size options.')
parser.add_argument('eds_file', type=argparse.FileType('r'), help='The EDS file to read')
parser.add_argument('num_patterns', type=int, help='How many patterns to search for')
parser.add_argument('--min-size', type=int, default=10, help='The minimum pattern size, default: 10')
parser.add_argument('--max-size', type=int, default=1000, help='The maximum pattern size, default: 1000')
parser.add_argument('--recreate-files', action='store_true', help='Delete and recreate the pattern file')
args = parser.parse_args()

if args.min_size > args.max_size or args.min_size <= 1:
	print 'Error: Invalid pattern lengths provided!'
	sys.exit(1)
if args.num_patterns <= 0:
	print 'Error: Invalid number of patterns!'
	sys.exit(1)

edsFileName = os.path.basename(args.eds_file.name)
edsFileName = ".".join(edsFileName.split('.')[0:-1])
patternsFile = './patterns/' + edsFileName + '_' + str(args.num_patterns)
if args.min_size == args.max_size:
	patternsFile += '_' + str(args.min_size)
else:
	patternsFile += '_' + str(args.min_size) + '-' + str(args.max_size)
patternsFile += '.txt'
resultsFile = './results/' + edsFileName + '_' + str(args.num_patterns)
if args.min_size == args.max_size:
	resultsFile += '_' + str(args.min_size)
else:
	resultsFile += '_' + str(args.min_size) + '-' + str(args.max_size)
resultsFile += '.txt'

if args.recreate_files:
	if os.path.isfile(patternsFile):
		os.remove(patternsFile)
	if os.path.isfile(resultsFile):
		os.remove(resultsFile)

#
# Create the patterns file
#
if not os.path.exists(patternsFile):
	maxPos = os.fstat(args.eds_file.fileno()).st_size - args.min_size - 1
	patterns = []
	while len(patterns) < args.num_patterns:
		args.eds_file.seek(0)
		args.eds_file.seek(random.randint(0, maxPos))
		rawPattern = args.eds_file.read(args.min_size + args.max_size)
		processedPattern = processRawPattern(rawPattern)
		if 'N' not in processedPattern and len(processedPattern) >= args.min_size:
			processedPattern = processedPattern[0:random.randint(args.min_size, args.max_size)]
			patterns.append(processedPattern)
			patterns = list(set(patterns))

	with open(patternsFile, 'w') as f:
		i = 0
		for p in patterns:
			if i != 0:
				f.write("\n")
			f.write(p)
			i += 1

#
# Search for the patterns with Multi-EDSM
#
cmd = './../../multiedsm -s ' + args.eds_file.name + ' -p ' + patternsFile
output = subprocess.Popen(['time ' + cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
with open(resultsFile, 'w') as f:
	f.write(output)
print output

#
# Check any patterns missing
#
patternIds = []
with open(resultsFile, 'r') as f:
	for line in f:
		if ',' in line and line[0] != 'P':
			patternIds.append(int(line.split(',')[-1]))

missingPatterns = set(range(args.num_patterns)) - set(patternIds)
if missingPatterns:
	with open(resultsFile, 'a') as f:
		f.write("\n\n")
		print "\n"
		for id in missingPatterns:
			missing = 'Failed to find pattern id ' + str(id)
			f.write(missing)
			print missing
