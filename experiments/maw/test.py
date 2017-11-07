
import os
import sys
import vcf
import argparse
import subprocess

#
# Parse command line arguments and check input
#
parser = argparse.ArgumentParser(description='Find the Minimum Absent Words (MAWs) for \
the reference Human Genome sequence then confirm these using 1000 Human Genomes data.')
parser.add_argument('chromosome', help='The chromosome number (1..22, X, Y)')
parser.add_argument('--min-size', type=int, default=3, help='The minimum size of a MAW (default = 3)')
parser.add_argument('--max-size', type=int, default=12, help='The maximum size of a MAW (default = 12)')
parser.add_argument('--memory', type=int, default=4096, help='How much RAM to use in MB, default=4096')
args = parser.parse_args()

if args.min_size > args.max_size:
	print 'Error: Invalid pattern sizes provided - min must be less than max!'
	sys.exit(1)
if args.min_size < 3:
	print 'Error: Invalid min size - must be 3 or greater!'
	sys.exit(1)
if args.chromosome not in ([str(x) for x in range(1, 23)] + ['X', 'Y']):
	print 'Error: Invalid chromosome provided - must be 1..22, X or Y'
	sys.exit(1)

humanGenome   = './references/hg37.fa'
emMAWOutput   = './em-maw-output/hg37.maws'
edsFile       = './eds/%s.eds' % args.chromosome
patternsFile  = './patterns/filtered_hg37.maws'
refFile       = './references/Homo_sapiens.GRCh37.75.dna.chromosome.%s.fa' % args.chromosome
resultsFile   = './results/%s_multiedsm.txt' % args.chromosome
falseMAWsFile = './results/falseMAWs_%s.txt' % args.chromosome
vcfgzFile = ''
if args.chromosome == 'Y':
	vcfgzFile = './vcf/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz'
elif args.chromosome == 'Y':
	vcfgzFile = './vcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
else:
	vcfgzFile = './vcf/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % args.chromosome

#
# Step 1: Run em-MAW on human genome
#

if not os.path.exists(emMAWOutput):
	cmd = './../../../maw/em-maw/em-maw -a DNA -i %s -o %s -k %d -K %d -m %d' % (humanGenome, emMAWOutput, args.min_size, args.max_size, args.memory)
	print subprocess.Popen([cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

#
# Step 2: Filter the MAWs file (remove MAWS containing Ns) to create a searchable patterns file
#

if not os.path.exists(emMAWOutput):
	print 'Error: could not locate MAW list.'
	sys.exit(1)
numMawsExtracted = 0
M = 0
if os.path.exists(patternsFile):
	with open(patternsFile, 'r') as pf:
		for line in pf:
			line = line.strip()
			numMawsExtracted += 1
			M += len(line)
else:
	with open(patternsFile, 'w') as pf:
		with open(emMAWOutput, 'r') as ef:
			ef.readline() # skip the > header
			for line in ef:
				line = line.strip()
				if not (len(line) == 0 or 'N' in line):
					if numMawsExtracted != 0:
						pf.write('\n')
					pf.write(line)
					numMawsExtracted += 1
					M += len(line)
print '%d MAWs of total length %d extracted for searching.' % (numMawsExtracted, M)

#
# Step 3: Search for the patterns with Multi-EDSM
#

output = ''
if numMawsExtracted == 0:
	print 'Not running Multi-EDSM on empty patterns file.'
	sys.exit(0)
if not os.path.exists(resultsFile):
	print 'Running Multi-EDSM...'
	cmd = './../../multiedsm -m %dm -s %s -p %s' % (args.memory, edsFile, patternsFile)
	output = subprocess.Popen(['time ' + cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
	with open(resultsFile, 'w') as f:
		f.write(output)
else:
	with open(resultsFile, 'r') as f:
		output = ''.join(f.readlines())
print output

#
# Step 4: extract patterns found
#

atMatches = False
positions = []
patternIds = []
patterns = {} # {patId:pattern,...}
with open(patternsFile, 'r') as f:
	for line in output.split('\n'):
		line = line.strip()
		if line == '':
			continue
		elif line.startswith('----'):
			atMatches = True
		elif atMatches and ',' in line:
			(pos, patId) = [int(x) for x in line.split(',')]
			positions.append(int(pos))
			patternIds.append(int(patId))
			if patId not in patterns:
				f.seek(0)
				i = 0
				while i < patId:
					f.readline()
					i += 1
				patterns[patId] = f.readline().strip()

print 'Verification of up to %d false positive MAWs starting.' % len(patterns.keys())

#
# Step 5: validate each match to see if it really exists and isn't a false positive
#

def processSample(POS, REF, ALTs, refPattern, sample, sampleStore, record, position, pattern):
	if sample.data.GT != None:
		for altIdx in sample.data.GT.split('|'):
			if altIdx == '0' or altIdx == '.' or '/' in altIdx:
				continue
			sampleALT = ALTs[int(altIdx)-1]
			if sampleALT.startswith('<') or sampleALT.startswith('.'):
				continue
			sampleName = sample.sample
			print record.POS, sampleALT, record.REF, position, refPattern, pattern
			break
			# if sampleName not in sampleStore:
			# 	sampleStore[sampleName] = [refPattern]
			# currMotifs = sampleStore[sampleName][:]
			# for motif in currMotifs:
			# 	temp = []
			# 	if POS < 1:
			# 		temp = list(REF)
			# 		superfluousPrefixLen = POS + len(REF)
			# 		for _ in range(superfluousPrefixLen):
			# 			try:
			# 				del temp[0]
			# 			except IndexError as e:
			# 				print _, POS, REF, sampleALT, superfluousPrefixLen, position, refPattern, pattern
			# 		temp = [sampleALT] + temp
			# 	else:
			# 		temp = list(motif)
			# 		try:
			# 			temp[POS] = sampleALT
			# 		except IndexError as e:
			# 			print POS, temp, motif, REF
			# 		for _ in range(0, min(len(REF)-1, len(refPattern)-POS)): # handle indels (reflen>1)
			# 			print 'Deleting at ' + position, len(REF)-1, len(refPattern), POS
			# 			del temp[POS+1]
			# 	temp = ''.join(temp)
			# 	sampleStore[sampleName].append(temp)

# confirmedPatternMatchIds = []
# vcf_reader = vcf.Reader(filename=vcfgzFile)
# for i in range(len(positions)):
# 	position = positions[i] + 1 # genome position starts at 1
# 	patternId = patternIds[i]
# 	pattern = patterns[patternId]
#
# 	# avoid repeat checking
# 	if patternId in confirmedPatternMatchIds:
# 		continue
#
# 	##
# 	# Open the reference genome and grab the ref pattern
# 	#
# 	refPattern = ''
# 	seq = ''
# 	with open(refFile, 'r') as rf:
# 		j = 0
# 		for line in rf:
# 			line = line.strip()
# 			if line == '' or line.startswith('>'):
# 				continue
# 			seq += line
# 			j += len(line)
# 			if j >= position:
# 				break
# 	refPattern = seq[position-len(pattern) : position]
# 	seq = ''
#
# 	##
# 	# Apply all the variants to the reference sequence at this position for any
# 	# individuals (sample) which have variants in this range
# 	#
# 	# for each variant in this range
# 	#	for each individual with an alternate allele
# 	#		apply the alt allele
# 	# Check if any of the sample's alleles match the MAW pattern.
# 	#
# 	vcfRecords = vcf_reader.fetch(args.chromosome, position - len(pattern) + 1, position)
# 	if vcfRecords:
# 		sampleStore = {}
# 		for record in vcfRecords:
# 			REF = str(record.REF)
# 			if REF[0] == '.' or REF[0] == '<':
# 				continue
# 			ALTs = [str(x) for x in record.ALT]
# 			POS = record.POS - (position - len(pattern) + 1)
# 			for sample in record.samples:
# 				processSample(POS, REF, ALTs, refPattern, sample, sampleStore, record, position, pattern)
# 				break
#
# 		matchFound = False
# 		for sampleId, alleles in sampleStore.items():
# 			for allele in alleles:
# 				if pattern in allele:
# 					matchFound = True
# 					confirmedPatternMatchIds.append(patternId)
# 					break
# 			if matchFound:
# 				break
#
# 	else:
# 		print 'Unexpected match for pattern #%d %s with %s' % (patternId, pattern, refPattern)
# 		if pattern == refPattern:
# 			confirmedPatternMatchIds.append(patternId)
















confirmedPatternMatchIds = []
vcf_reader = vcf.Reader(filename=vcfgzFile)
for i in range(len(positions)):
	matchEndPosition = positions[i] + 1 # genome position starts at 1
	patternId = patternIds[i]
	mawPattern = patterns[patternId]
	matchStartPosition = matchEndPosition - len(mawPattern) + 1

	# avoid repeat checking
	if patternId in confirmedPatternMatchIds:
		continue

	##
	# Open the reference genome and grab the ref pattern
	#
	refPattern = ''
	seq = ''
	with open(refFile, 'r') as rf:
		j = 0
		for line in rf:
			line = line.strip()
			if line == '' or line.startswith('>'):
				continue
			seq += line
			j += len(line)
			if j >= matchEndPosition:
				break
	refPattern = seq[matchEndPosition-len(mawPattern) : matchEndPosition]
	seq = ''

	##
	# Apply all the variants to the reference sequence at this position for any
	# individuals (sample) which have variants in this range
	#
	sampleStore = {}
	for record in vcf_reader.fetch(args.chromosome, matchStartPosition, matchEndPosition):
		POSstart = record.POS
		REF = str(record.REF)
		REFlen = len(REF)
		POSend = POSstart + REFlen - 1
		if POSend != 47771698:
			continue
		if REF[0] == '.' or REF[0] == '<':
			continue
		ALTs = [str(x) for x in record.ALT]
		for sample in record.samples:
			if sample.data.GT != None:
				for altIdx in sample.data.GT.split('|'):
					if '/' in altIdx or altIdx == '0' or altIdx == '.':
						continue
					sampleALT = ALTs[int(altIdx)-1]
					if sampleALT.startswith('<') or sampleALT.startswith('.'):
						continue
					sampleName = sample.sample
					if sampleName not in sampleStore:
						sampleStore[sampleName] = [refPattern]
					currMotifs = sampleStore[sampleName][:]
					indexToSubstitute = POSstart - matchStartPosition

					for motif in currMotifs:
						temp = list(motif)
						if indexToSubstitute < 0:
							delTimes = POSend - matchStartPosition
							print 'underflow:', delTimes, POSend, matchStartPosition, motif, sampleALT, mawPattern, refPattern
							for _ in range(min(delTimes, len(motif) - 1)):
								del temp[0]
							temp = [sampleALT] + temp
							print temp
						elif REFlen > 1 and POSend > matchEndPosition:
							delTimes = matchEndPosition - POSstart
							print 'overflow:', delTimes, POSend, matchStartPosition, motif, sampleALT, mawPattern, refPattern, REFlen
							for _ in range(delTimes - 1):
								del temp[indexToSubstitute]
							temp = temp + [sampleALT]
						else:
							print 'simpleSub:', mawPattern, indexToSubstitute, sampleALT, refPattern #, temp
							temp[indexToSubstitute] = sampleALT
							for _ in range(REFlen - 1):
								del temp[indexToSubstitute + 1]
						temp = ''.join(temp)
						sampleStore[sampleName].append(temp)

	matchFound = False
	for sampleId, alleles in sampleStore.items():
		for allele in alleles:
			if mawPattern in allele:
				matchFound = True
				confirmedPatternMatchIds.append(patternId)
				break
		if matchFound:
			break

#
# Step 6: Output the confirmedPatternMatchIds (false MAWs) list - this will be
# used later to reduce the size of the original MAWs list -- see summary.py
#

if len(confirmedPatternMatchIds) == 0:
	print 'Verified all MAWs are valid and matches found by MultiEDSM belong to multiple independent samples'
with open(falseMAWsFile, 'w') as f:
	i = 0
	for patternId in confirmedPatternMatchIds:
		print 'False MAW #%s (%s) identified.' % (str(patternId), patterns[patternId])
		if i != 0:
			f.write('\n')
		f.write(str(patternId))
		i += 1
print 'All done!'
