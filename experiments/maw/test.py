
import os
import sys
import argparse
import subprocess

#
# Parse command line arguments and check input
#
parser = argparse.ArgumentParser(description='Find the Minimum Absent Words (MAWs) \
in a given reference sequence and MAW length and output these motifs to a file. \
Then remove invalid entries from the file to create a patterns file. Search the \
given Pan-Genome (EDS file) and note all matches found for verification to see \
if it is really a MAW. Note: Output files will be named after the EDS file and \
put in the results folder and named {true|false}MAWs_{chr}_{pattern_size}.txt')
parser.add_argument('ref_file', type=argparse.FileType('r'), help='The Reference Chromosome file to look for MAWs in')
parser.add_argument('eds_file', type=argparse.FileType('r'), help='The EDS Pan-Genome file to search in')
parser.add_argument('pattern_size', type=int, help='The motif size of MAWs to discover')
parser.add_argument('--recreate-files', action='store_true', help='Delete and recreate the pattern file')
parser.add_argument('--memory', type=int, default=4096, help='How much RAM to use in MB, default=4096')
args = parser.parse_args()

if args.pattern_size < 3 or args.pattern_size > 13:
	print 'Error: Invalid pattern lengths provided! Must be between 3 and 13 inclusive!'
	sys.exit(1)

refFileName = os.path.basename(args.ref_file.name)
edsFileName = os.path.basename(args.eds_file.name)
edsFileName = ".".join(edsFileName.split('.')[0:-1])
emMAWOutput = './em-maw-output/' + edsFileName + '_' + str(args.pattern_size) + '.txt'
patternsFile = './patterns/' + edsFileName + '_' + str(args.pattern_size) + '.txt'
resultsFile = './results/' + edsFileName + '_' + str(args.pattern_size) + '.txt'
if edsFileName == 'Y':
	vcfFile = './vcf/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz'
else:
	vcfFile = './vcf/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % edsFileName
verificationFile = './results/verification_' + edsFileName + '_' + str(args.pattern_size) + '.txt'
trueMAWsFile = './results/trueMAWS_' + edsFileName + '_' + str(args.pattern_size) + '.txt'
falseMAWsFile = './results/falseMAWS_' + edsFileName + '_' + str(args.pattern_size) + '.txt'

try:
	subprocess.Popen(['vcftools --help'], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()
except OSError as e:
	if e.errno == os.errno.ENOENT:
		print 'Error: It appears you do not have VCFTools installed!'
	else:
		print e
	sys.exit(1)

if args.recreate_files:
	if os.path.isfile(emMAWOutput):
		os.remove(emMAWOutput)
	if os.path.isfile(patternsFile):
		os.remove(patternsFile)
	if os.path.isfile(resultsFile):
		os.remove(resultsFile)
	if os.path.isfile(verificationFile):
		os.remove(verificationFile)
	if os.path.isfile(trueMAWsFile):
		os.remove(trueMAWsFile)
	if os.path.isfile(falseMAWsFile):
		os.remove(falseMAWsFile)
	for fileName in os.listdir('./references'):
		if fileName.startswith(refFileName) and (fileName.endswith('.bwt5') or fileName.endswith('.lcp5') or fileName.endswith('.sa5')):
			os.remove('./references/' + fileName)
	for fileName in os.listdir('./vcf'):
		if fileName.endswith('.recode.vcf'):
			os.remove('./vcf/' + fileName)

#
# Step 1: Search for MAWs in the reference sequence
#
if not os.path.exists(emMAWOutput):
	cmd = './../../../maw/em-maw/em-maw -a DNA -i ' + args.ref_file.name + ' -o ' \
		+ emMAWOutput + ' -k ' + str(args.pattern_size) + ' -K ' + str(args.pattern_size) \
		+ ' -m ' + str(args.memory)
	for fileName in os.listdir('./references'):
		if refFileName + '_' in fileName:
			cmd += ' -c 0'
			break
	output = subprocess.Popen([cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
	print output

#
# Step 2: Filter the MAWs file to create a patterns file
#
if not os.path.exists(emMAWOutput):
	print 'Error: could not locate MAW list.'
	sys.exit(1)
numMawsExtracted = 0
patterns = []
with open(emMAWOutput, 'r') as f:
	f.readline() # skip the > header
	for line in f:
		line = line.strip()
		if not (len(line) == 0 or ('N' in line)):
			patterns.append(line)
			numMawsExtracted += 1
if numMawsExtracted > 0:
	with open(patternsFile, 'w') as pf:
		pf.write('\n'.join(patterns))
print '%d MAWs extracted for searching.' % numMawsExtracted

#
# Step 3: Search for the patterns with Multi-EDSM
#
output = ''
if numMawsExtracted == 0:
	print 'Not running Multi-EDSM on empty patterns file.'
	sys.exit(0)
elif not os.path.exists(resultsFile):
	print 'Running Multi-EDSM...'
	cmd = './../../multiedsm -m ' + str(args.memory) + 'm -s ' + args.eds_file.name + ' -p ' + patternsFile
	output = subprocess.Popen(['time ' + cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
	with open(resultsFile, 'w') as f:
		f.write(output)
	print output
else:
	with open(resultsFile, 'r') as f:
		output = '\n'.join(f.readlines())
print '\nMultiEDSM search completed. Now creating the VCF files for the matches.'

#
# Step 4: Read the positions from the MultiEDSM output and for each position
# matched create a VCF file containing the mutations in the pattern with VCFTools
#
vcfFiles = []
matchedPatterns = []
atMatches = False
for line in output.split('\n'):
	if line == '':
		continue
	elif line.startswith('----'):
		atMatches = True
	elif atMatches and ',' in line:
		(position, patternId) = line.split(',')
		position = int(position) + 1 #genomes start at position 1
		patternId = int(patternId)
		vcfPosFileName = './vcf/' + edsFileName + '_' + str(position) + '_' + str(patternId) + '.recode.vcf'
		if not os.path.isfile(vcfPosFileName):
			cmd = 'vcftools --chr ' + edsFileName + ' --from-bp ' + str(position - args.pattern_size) \
			    + ' --to-bp ' + str(position) + ' --gzvcf ' + vcfFile + ' --out ./vcf/' \
			    + edsFileName + '_' + str(position) + '_' + str(patternId) + ' --recode'
			subprocess.Popen([cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()
		vcfFiles.append(vcfPosFileName)
		matchedPatterns.append({'pos':position, 'pat':patterns[patternId]})
if not atMatches:
	print 'There were no matches found with Multi-EDSM. Quitting.'
	sys.exit(0)
else:
	print 'Finished creating the VCF files for each MAW discovered in the Pan-Genome'

#
# Step 5: Read the newly created VCF files of each match and extract the sample
# ids and associate mutations with them at the correct position. Verify sample
# DNA matches the MAW indicating it's not a true MAW
#

trueMAWs = []
falseMAWs = []
falseMAWSamples = []
for id in range(len(matchedPatterns)):
	matchPos = matchedPatterns[id]['pos'] # 0-based index
	matchPat = matchedPatterns[id]['pat'] # MAW
	vcfFile = vcfFiles[id]

	maw = matchPat
	pos = matchPos
	ref = ''
	maw_size = len(maw)

	##
	# Open the reference genome and grab the ref pattern
	#
	seq = ''
	i = 0
	args.ref_file.seek(0)
	for line in args.ref_file:
		line = line.strip()
		if line == '' or line.startswith('>'):
			continue
		seq += line
		i += len(line)
		if i > pos:
			break
	ref = seq[pos-maw_size : pos]
	seq = ''

	##
	# Open up the VCF file and grab the mutations matching the MAW
	#
	sampleIds = []
	samples = {}
	with open(vcfFile, 'r') as vcf_file:
		for line in vcf_file:
			if line.startswith('##'):
				continue
			elif line.startswith('#'):
				sampleIds = [x.strip() for x in line.split('\t')]
				sampleIds = sampleIds[9:] #remove the CHROM, POS, ID... cols leaving just the sample IDs
			else:
				record = [x.strip() for x in line.split('\t')]
				(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) = record[0:9]
				POS = int(POS)
				if POS < (pos - maw_size) or POS > pos:
					print 'Warning: For MAW ' + maw + ' skipped ' + ID + ' at position: ' + str(POS)
					continue
				if '.' in REF + ALT or '<' in REF + ALT:
					continue
				record = record[9:] #remove the CHROM, POS, ID... cols leaving just the sample IDs
				POS = POS - (pos - maw_size)
				reflen = len(REF)
				if maw[POS-1 : POS-1+reflen] == ref[POS-1 : POS-1+reflen]:
					continue
				altAlleles = ALT.split(',')
				for i, v in enumerate(record):
					for altIndex in v.split('|'):
						if altIndex == '0' or altIndex == '.':
							continue
						alt = altAlleles[int(altIndex) - 1]
						if alt == maw[POS-1 : POS-1+len(alt)]:
							if POS not in samples:
								samples[POS] = [sampleIds[i]]
							else:
								samples[POS].append(sampleIds[i])

	##
	# Confirm one or more sample ids have all the mutations of the MAW
	#
	sampleSets = []
	guiltySamples = []
	for k, v in samples.items():
		sampleSets.append(set(v))
	if len(sampleSets) > 0:
		guiltySamples = set.intersection(*sampleSets)
	if len(guiltySamples) == 0:
		print 'MAW ' + maw + ' is valid'
		trueMAWs.append(maw)
	else:
		print 'False MAW ' + maw + ' found at position ' + str(pos + 1) + ' was ' \
			+ 'discovered for samples: ' + ', '.join(guiltySamples)
		falseMAWs.append(maw)
		falseMAWSamples.append(list(guiltySamples))
print 'Finished verifying the MAWs discovered!'

#
# Step 6: Output the true and false MAWs
#
if len(trueMAWs) > 0:
	with open(trueMAWsFile, 'w') as f:
		for i, maw in enumerate(trueMAWs):
			if i != 0:
				f.write('\n')
			f.write(maw)
else:
	print 'No true MAWS in discovered!'
if len(falseMAWs) > 0:
	with open(falseMAWsFile, 'w') as f:
		for i, maw in enumerate(falseMAWs):
			if i != 0:
				f.write('\n')
			f.write(maw)
else:
	print 'No false MAWs discovered!'

print ''
print str(len(trueMAWs)) + ' potentially true MAWs confirmed out of ' + str(numMawsExtracted) + ':'
print '\n'.join(trueMAWs)
