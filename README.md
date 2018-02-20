## Multi-EDSM

Multiple string matching in Elastic Degenerate Strings.

If you make use of this software, please cite the following:

> S. P. Pissis, A. Retha, "Dictionary Matching in Elastic-Degenerate Texts with Applications in Searching VCF Files On-line", (submitted to SEA2018).

### How to Install

You may need to install one or more libraries before compiling Multi-EDSM. Please read INSTALL.md.

### How to Execute

If you are searching an EDS format file:

`~$ ./multiedsm --sequence-file seq.eds --patterns-file patterns.txt --mem-limit 4g`

If you are searching a reference sequence in FASTA format along with its VCF file:

`~$ ./multiedsm --sequence-file reference.fasta --variants-file variants.vcf --patterns-file patterns.txt --mem-limit 4g`

For more options read the help instructions: `~$ ./multiedsm --help`.

The patterns.txt file should be a set of strings separated by line breaks.

If you want to use a compressed vcf file (*.vcf.gz), please make sure its accompanying tbi file is also present in the same directory. You can also use `Tabix` to generate a tbi file.

### License

GNU GPLv3 License; Copyright (C) 2018 Solon P. Pissis and Ahmad Retha.
