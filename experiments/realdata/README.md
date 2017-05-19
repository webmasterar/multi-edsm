## Obtaining realdata experiment files and instructions

The realdata experiment files are from the 1000 Genomes Project. You will need to download three sets of files into the datasets folder to run the experiment.

The reference sequences were downloaded from: ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/

For example, we download the sequence file for chromosome Y:

`~$ wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.Y.fa.gz`

The reference sequence files must be extracted to work with EDSM:

`~$ gzip -d Homo_sapiens.GRCh37.75.dna.chromosome.Y.fa.gz`

The variant call files (vcf.gz) and their associated tbi files were downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

Download the vcf and tbi files for chromosome Y:

`~$ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz`

`~$ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi`

Then you can use the EDSO tool to convert the reference sequence fasta file with its VCF file to EDS format. EDSO can be obtained from: `https://github.com/webmasterar/edso`.

With the EDS file created, you will be able to run the test.py script.
