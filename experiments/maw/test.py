
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
		output = '\n'.join(f.readlines())
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

confirmedPatternMatchIds = []
vcf_reader = vcf.Reader(filename=vcfgzFile)
for i in range(len(positions)):
	position = positions[i] + 1 # genome position starts at 1
	patternId = patternIds[i]
	pattern = patterns[patternId]

	# avoid repeat checking
	if patternId in confirmedPatternMatchIds:
		continue

	##
	# Open the reference genome and grab the ref pattern
	#
	ref = ''
	seq = ''
	with open(refFile, 'r') as rf:
		j = 0
		for line in rf:
			line = line.strip()
			if line == '' or line.startswith('>'):
				continue
			seq += line
			j += len(line)
			if j >= position:
				break
	ref = seq[position-len(pattern) : position]
	seq = ''

	##
	# Apply all the variants to the reference sequence at this position for any
	# individuals (samples) which have variants in this range
	#
	# for each variant in this range
	#	for each individual with an alternate allele
	#		apply the alt allele
	#
	vcfRecords = vcf_reader.fetch(args.chromosome, position - len(pattern), position)
	if vcfRecords:
		samples = {}
		for record in vcfRecords:
			if record.REF[0] == '.' or record.REF[0] == '<':
				continue
			if '.' in record.ALT or '<' in record.ALT:
				continue
			POS = record.POS - (position - len(pattern))
			for sample in record.samples:
				if sample.data.GT == None:
					continue
				sampleName = sample.sample
				for altIndex in sample.data.GT.split('|'):
					if altIndex == '0' or altIndex == '.' or '/' in altIndex:
						continue
					if sampleName not in samples:
						samples[sampleName] = [ref]
					sampleALT = str(record.ALT[int(altIndex) - 1])
					if sampleALT.startswith('<'):
						continue
					currMotifs = samples[sampleName][:]
					for motif in currMotifs:
						temp = list(motif)
						temp[POS-1] = sampleALT
						temp = ''.join(temp)
						samples[sampleName].append(temp)

		##
		# Check if any of the sample's alleles match the MAW pattern
		#
		matchFound = False
		for sampleId, alleles in samples.items():
			for allele in alleles:
				if pattern in allele:
					matchFound = True
					confirmedPatternMatchIds.append(patternId)
					break
			if matchFound:
				break

	else:
		print 'Unexpected match for pattern %d %s with %s' % (patternId, pattern, ref)
		if pattern == ref:
			confirmedPatternMatchIds.append(patternId)


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
