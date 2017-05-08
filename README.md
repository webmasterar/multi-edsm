## Multi-EDSM

Multiple string matching in Elastic Degenerate Strings.

If you make use of this software, please cite the following:

> Robert Grossi, Costas S. Iliopoulos, Chang Liu, Nadia Pisanti, Solon P. Pissis, Ahmad Retha, Giovanna Rosone, Fatima Vayani and Luca Versari; On-line pattern matching on similar texts; (in preparation).

### How to Install

You may need to install one or more libraries before compiling Multi-EDSM. Please read INSTALL.md.

### How to Execute

`~$ ./multiedsm seq.txt patterns.txt`

or

`~$ ./multiedsm reference.fasta variants.vcf pattern.txt`

The patterns.txt file should be a set of strings separated by line breaks.

If you want to use a compressed vcf file (*.vcf.gz), please make sure its accompanying tbi file is also present in the same directory. You can also use `Tabix` to generate a tbi file.

### License

GNU GPLv3 License; Copyright (C) 2017 Chang Liu, Solon P. Pissis, Ahmad Retha and Fatima Vayani.
