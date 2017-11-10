## 1000 Human Genomes Project MAW Extraction & Verification

### Software required

To run the experiment requires you obtain some additional software:-

First, clone the Github project [MAW](https://github.com/solonas13/maw) to
./../../../maw then go into the *em-maw* directory and read the instructions on
how to build it.

Secondly, we use [PyVCF](https://github.com/jamescasbon/PyVCF) python module to
read the VCF files. Go ahead and install that. This is usually done by running
`pip install pyvcf`. The experiment was run under Python v2.7 using PyVCF 0.6.7.

Thirdly, we create *EDS* (Elastic Degenerate String) files from the FASTA+VCF
files using the Github project [EDSO](https://github.com/webmasterar/edso).
Clone it to ./../../../edso then read the installation instructions to build it.

Finally, confirm that you have built the *Multi-EDSM* project and the binary
./../../multiedsm exists.

### Obtaining files

First, download the data files (FASTA+VCF) you can execute
`./download1000GenomesFiles.sh`. This downloads about 18GB of files.

Then you need to create the combined Human reference genome by running
`./createCombinedGenome.sh`. This will require 6.3GB additional disk space.

Then you have to create the EDS files. This is a very slow process taking days
and requires an additional 3.5GB of disk space. Execute `./createEDSFiles.sh`.

### Running the experiment

To start the experiment simply run: `./run-experiment.sh`. This will run *test.py*
for each of the chromosomes. To see the options available try running
`python test.py --help`.

Executing *run-experiment.sh* is a very slow process and will take days to complete.
The first command executed in the file will run *em-maw* to create a list of MAWs
identified from the combined reference Human genome and will require 34GB of disk
space to create an index. Subsequent commands will skip this initial indexing process.
Each command will run *Multi-EDSM* to find any false-positive MAWs in each
chromosome and verify they correspond to an individual (sample) from the 1000
Human Genomes project.

### Results

The `./results` folder contains the experiment results. Files named
*%_multiedsm.txt* contain the results of running Multi-EDSM and contain potential
false MAWs. Files named *falseMAWs_%.txt* contain confirmed false MAWs.

A full list of the verified MAWs is obtained by running `python summary.py`.
The confirmed MAWs are listed in `./results/confirmedMAWs.maws`.
