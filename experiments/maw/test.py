
import os
import sys
import argparse
import subprocess

#
# Parse command line arguments and check input
#
parser = argparse.ArgumentParser(description='Find the Minimum Absent Words (MAWs) \
in a given reference sequence and MAW length and output these motifs to a file. Then \
remove invalid entries from the file to create a patterns file. Search the given \
Pan-Genome (EDS file) and note all matches found for verification to see if it is \
really a MAW. Note: Output files will be named after the EDS file.')
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
	for fileName in os.listdir('./references'):
		if fileName.startswith(refFileName) and (fileName.endswith('.bwt5') or fileName.endswith('.lcp5') or fileName.endswith('.sa5')):
			os.remove('./references/' + fileName)

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
i = 0
patterns = []
with open(emMAWOutput, 'r') as f:
	f.readline() # skip the > header
	for line in f:
		line = line.strip()
		if not (len(line) == 0 or ('N' in line)):
			patterns.append(line)
			i += 1
if i > 0:
	with open(patternsFile, 'w') as pf:
		pf.write('\n'.join(patterns))
print '%d MAWs extracted for searching.' % i

#
# Step 3: Search for the patterns with Multi-EDSM
#
output = ''
if i == 0:
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

#
# Step 4: Read the positions from the MultiEDSM output and for each position
# matched create a VCF file containing the mutations in the pattern with VCFTools
#
vcfFiles = []
matchPatterns = []
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
		cmd = 'vcftools --chr ' + edsFileName + ' --from-bp ' + str(position - args.pattern_size) \
		    + ' --to-bp ' + str(position) + ' --gzvcf ' + vcfFile + ' --out ./vcf/' \
		    + edsFileName + '_' + str(position) + '_' + str(patternId) + ' --recode'
		subprocess.Popen([cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()
		vcfFiles.append('./vcf/' + edsFileName + '_' + str(position) + '_' + str(patternId) + '.recode.vcf')
		matchPatterns.append({'pos':position, 'pat':patterns[patternId]})
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
def getRefPatternEndingAtPos(matchPosition, pattern_size):
	args.ref_file.seek(0)
	seq = ''
	i = 0
	for line in args.ref_file:
		line = line.strip()
		if line == '' or line.startswith('>'):
			continue
		seq += line
		i += len(line)
		if i > matchPosition:
			break
	return seq[matchPosition - pattern_size: matchPosition]

for vcfId in range(len(vcfFiles)):
	with open(vcfFiles[vcfId], 'r') as v:
		#
		# 5a: Applying SNPs to each sample
		#
		matchPattern = matchPatterns[vcfId]['pat'] # MAW
		matchPosition = matchPatterns[vcfId]['pos']
		refPattern = getRefPatternEndingAtPos(matchPosition, args.pattern_size) # pattern in ref sequence
		sampleIds = []
		samples = {}
		atRecord = False
		for line in v:
			if line.startswith('##'):
				continue
			elif line.startswith('#'):
				sampleIds = [x.strip() for x in line.split('\t')]
				sampleIds = sampleIds[9:] #remove the CHROM, POS, ID... cols leaving just the sample IDs
				atRecord = True
			elif atRecord:
				record = [x.strip() for x in line.split('\t')]
				(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) = record[0:9]
				POS = int(POS) - (matchPosition - args.pattern_size)
				if '.' in REF + ALT or '<' in REF + ALT:
					continue
				print str(POS), str(matchPosition)
				record = record[9:] #remove the CHROM, POS, ID... cols leaving just the sample alleles
				#for colId in range(len(record)):
				#	if record[colId] != '0':
				#		key = sampleIds[colId]
				# 		if key not in samples:
				#			samples[key] = refPattern
				# 		temp = list(samples[key])
				#		temp[POS - 1] = ALT
				#		samples[key] = ''.join(temp)
				#		print ALT + ' ' + key + ': ' + samples[key] + ' <= MAW:' + matchPattern + ' ref:' + refPattern
				# OR:
				# for colId in range(len(record)):
				# 	if not (record[colId] == '0' or record[colId] == '0|0'):
				# 		key = sampleIds[colId]
				# 		altAlleles = ALT.split(',')
				# 		altAlleleIndices = [int(x) for x in record[colId].split('|')]
				# 		if key not in samples:
				# 			samples[key] = [refPattern]
				# 		for altIdx in altAlleleIndices:
				# 			if altIdx != 0:
				# 				altIdx -= 1
				# 				temp = list(samples[key][0])
				# 				temp[POS - 1] = altAlleles[altIdx]
				# 				temp = ''.join(temp)
				# 				samples[key].append(temp)
				# 				print altAlleles[altIdx] + ' ' + key + ': ' + temp + ' <= MAW:' + matchPattern + ' ref:' + refPattern
				# OR:
				for colId in range(len(record)):
					if not (record[colId] == '0' or record[colId] == '0|0' or record[colId] == '.'):
						key = sampleIds[colId]
						altAlleles = ALT.split(',')
						altAlleleIndices = [int(x) for x in record[colId].split('|')]
						if key not in samples:
							samples[key] = [refPattern]
						currMotifs = samples[key][:]
						for motif in currMotifs:
							for altIdx in altAlleleIndices:
								if altIdx != 0:
									altIdx -= 1
									temp = list(motif)
									temp[POS - 1] = altAlleles[altIdx]
									temp = ''.join(temp)
									samples[key].append(temp)
									print altAlleles[altIdx] + ' ' + key + ': ' + temp + ' <= MAW:' + matchPattern + ' ref:' + refPattern

		#
		# 5b: Verifying MAW pattern doesn't match any sample
		#
		#falseMAWcount = 0
		#for key, val in samples.items():
		#	if matchPattern in val:
		#		falseMAWcount += 1
		#		msg = val + ' false MAW found for sample ' + key + ' at position ' + str(matchPosition)
		#		print msg
		#		with open(verificationFile, 'a') as v:
		#			v.write(msg + '\n')
		#print '%d false MAWs identified.' % falseMAWcount
		# OR:
		falseMAWcount = 0
		for key, motifArr in samples.items():
			for motif in set(motifArr):
				if matchPattern in motif:
					falseMAWcount += 1
					msg = motif + ' false MAW found for sample ' + key + ' at position ' + str(matchPosition)
					print msg
					with open(verificationFile, 'a') as v:
						v.write(msg + '\n')
		print '%d false MAWs identified for position:%d, pattern:%s.' % (falseMAWcount, matchPatterns[vcfId]['pos'], matchPatterns[vcfId]['pat'])
print 'MAW match verification complete.'
print 'All done!'
