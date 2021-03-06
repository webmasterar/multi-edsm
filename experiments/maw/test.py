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
parser.add_argument('--memory', type=int, default=4096, help='How much RAM to use in MB (default = 4096)')
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
elif args.chromosome == 'X':
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
 			M += len(line.strip())
 			numMawsExtracted += 1
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
#patterns = {} # {patId:pattern,...}
patterns = []
with open(patternsFile, 'r') as f:
	for line in f:
		patterns.append(line.rstrip())
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

print 'Verification of up to %d potential false positive MAWs starting.' % len(set(patternIds))

#
# Step 5: validate each match to see if the MAW really exists and therefore isn't
# really a MAW
#

refText = ''
with open(refFile, 'r') as rf:
	for line in rf:
		line = line.strip()
		if line == '' or line.startswith('>'):
			continue
		refText += line

##
# Get range of bases from reference fasta file. e.g. chr21:9411239-9411240
# @param begin 1-based index in reference sequence e.g. 9411239
# @param end 1-based index in reference up-to-but-not-including e.g. 9411240
# @return string 'G'
def getRefText(begin, end):
	if begin > end:
		begin, end = end, begin
	elif begin == end or begin < 1 or end < 1:
		return ''
	return refText[begin - 1 : end - 1]

disqualifiedMAWIds = []
vcf_reader = vcf.Reader(filename=vcfgzFile)
for i in range(len(positions)):
	matchEndPosition = positions[i] + 1 # genome position starts at 1 (1-based index)
	patternId = patternIds[i]
	mawPattern = patterns[patternId]
	matchStartPosition = matchEndPosition - len(mawPattern) + 1

	# avoid repeat checking
	if patternId in disqualifiedMAWIds:
		continue

	##
	# Grab the variants (alleles) in this range (matchStartPosition - matchEndPosition)
	# from the VCF file
	#
	requiredAdditionalPrefix = 0
	reffetchbegin = matchStartPosition
	reffetchend = matchEndPosition + 1
	alleles = [] # [{ref:'G', startPos:0, alt:['A','C'], samples:{'HG001': ['0','1'], 'HG002': ['0','0'], ...}]
	for record in vcf_reader.fetch(args.chromosome, matchStartPosition, matchEndPosition):
		POSstart = record.POS
		REF = str(record.REF)
		REFlen = len(REF)
		POSend = POSstart + REFlen - 1
		if REF[0] == '.' or REF[0] == '<':
			continue
		ALTs = []
		maxDel = 0
		for alt in record.ALT:
			alt = str(alt)
			if alt.startswith('<') or alt.startswith('.'):
				continue
			ALTs.append(alt)
			if len(alt) < REFlen and REFlen-len(alt) > maxDel:
				maxDel = REFlen-len(alt)
		if len(ALTs) == 0:
			continue
		requiredAdditionalPrefix += maxDel
		samples = {}
		for sample in record.samples:
			if sample.data.GT != None:
				samples[sample.sample] = [str(x) for x in str(sample.data.GT).replace('/', '|').split('|')] #even if not phased, treat it as if phased because subsequent records in the block are, see: https://www.biostars.org/p/109632/
		alleles.append({'ref':REF, 'startPos':POSstart, 'alt':ALTs, 'samples':samples})
		if POSstart < reffetchbegin:
			reffetchbegin = POSstart

	##
	# Open the reference chromosome and grab the ref sequence for the range
	#
	reffetchbegin = reffetchbegin - requiredAdditionalPrefix
	refSequence = getRefText(reffetchbegin, reffetchend)

	##
	# Recreate the sequence for every sample (person) that has an alternate allele
	# in this range and do so conserving chromosomal haplotype
	#
	phasedSampleStore = {} # {'HG001': ['GTCCA','GTACA'], 'HG002': ['GTCCA', 'GTCCA'], ...}
	sampleShiftIdx = {} # {'HG001': [0,0], 'HG002': [0,1], ...}
	for allele in alleles:
		if allele['ref'].startswith('.') or allele['ref'].startswith('<'):
			continue
		for sampleId, genotypes in allele['samples'].items():
			for haplotypeIdx, genotypeIdx in enumerate(genotypes):
				if genotypeIdx == '.':
					continue
				genotypeIdx = int(genotypeIdx)
				if genotypeIdx == 0:
					variant = allele['ref']
				else:
					variant = allele['alt'][genotypeIdx - 1]
					if variant.startswith('.') or variant.startswith('<'):
						continue
				if sampleId not in sampleShiftIdx:
					sampleShiftIdx[sampleId] = [0, 0]
				if sampleId not in phasedSampleStore:
					phasedSampleStore[sampleId] = [refSequence, refSequence]
				indexToSubstitute = allele['startPos'] - reffetchbegin + sampleShiftIdx[sampleId][haplotypeIdx]
				temp = phasedSampleStore[sampleId][haplotypeIdx]
				if indexToSubstitute >= len(temp) or indexToSubstitute < 0:
					continue
				temp = list(temp)
				temp[indexToSubstitute] = variant
				for j in range(1, min(len(allele['ref']), len(temp) - indexToSubstitute)):
					del temp[indexToSubstitute+1]
				phasedSampleStore[sampleId][haplotypeIdx] = ''.join(temp)
				sampleShiftIdx[sampleId][haplotypeIdx] += len(variant) - len(allele['ref'])

	##
	# Check if the recreated sequences of any sample's haplotype match the MAW
	# and add it to the disqualified list if it does
	#
	for sampleId, haplotypes in phasedSampleStore.items():
		if mawPattern in haplotypes[0] or mawPattern in haplotypes[1]:
			disqualifiedMAWIds.append(patternId)
			break

#
# Step 6: Output the disqualifiedMAWIds (false MAWs) list - this will be
# used later to reduce the size of the original MAWs list -- see summary.py
#

if len(disqualifiedMAWIds) == 0:
	print 'Verified all MAWs are valid and matches found by MultiEDSM belong to multiple independent samples.'
with open(falseMAWsFile, 'w') as f:
	i = 0
	for patternId in disqualifiedMAWIds:
		print 'False MAW #%s (%s) identified.' % (str(patternId), patterns[patternId])
		if i != 0:
			f.write('\n')
		f.write(str(patternId))
		i += 1
print 'All done!'
