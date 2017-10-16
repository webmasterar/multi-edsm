
import argparse
import subprocess

def main():
	parser = argparse.ArgumentParser(description='Given a MAW string and it\'s \
expected postion along with the associated VCF and Reference genome files, this \
program verifies the MAW is present (i.e. it\'s not really a MAW) and returns \
the sampleIds of the genomes the MAW exists in.')
	parser.add_argument('maw', help='The MAW to be verified')
	parser.add_argument('position', type=int, help='The position where the MAW was found (0-based index)')
	parser.add_argument('ref_file', type=argparse.FileType('r'), help='The Reference Chromosome file to look for MAWs in')
	parser.add_argument('vcf_file', type=argparse.FileType('r'), help='The VCF file to search in')
	args = parser.parse_args()

	maw = args.maw
	pos = args.position
	ref = ''
	maw_size = len(maw)

	##
	# Open the reference genome and grab the ref
	#
	seq = ''
	i = 0
	for line in args.ref_file:
		line = line.strip()
		if line == '' or line.startswith('>'):
			continue
		seq += line
		i += len(line)
		if i > pos:
			break
	ref = seq[pos-maw_size : pos]
	seq = None

	##
	# Open up the VCF file and grab the mutations matching the MAW
	#
	sampleIds = []
	samples = {}
	for line in args.vcf_file:
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
				print 'skipped ' + ID + ' at position: ' + str(POS)
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
	for k, v in samples.items():
		sampleSets.append(set(v))
	guiltySamples = set.intersection(*sampleSets)
	if len(guiltySamples) == 0:
		print 'MAW ' + maw + ' is valid'
	else:
		print 'False MAW ' + maw + ' found at position ' + str(pos + 1) + ' was ' \
			+ 'discovered for samples: ' + ', '.join(guiltySamples)

if __name__ == '__main__':
	main()
