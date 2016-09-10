# EMContaminator
Infer the fraction of sequencing data that comes from each individual in a population

#### Author: Thomas Willems <tfwillems@gmail.com> <br> License: GNU v3

## Installation
EMContaminator requires a standard c++ compiler and [CMake](http://www.cmake.org/download/).
To obtain it and all of its associated submodules, use:

    % git clone --recursive https://github.com/tfwillems/EMContaminator.git

To build, use Make:

    % cd EMContaminator
    % make

## Quick Start
To run EMContaminator, use the syntax:
```
./EMContaminator --bam <reads.bam> --vcf <snps.vcf.gz>
```

* **--bam**:   <reads.bam>              BAM file containing reads from sequencing run
* **--vcf**:  <snps.vcf.gz>            Bgzipped VCF file containing SNP genotypes for a population of individuals
