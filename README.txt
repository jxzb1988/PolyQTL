# PolyQTL: A Bayesian method to detect multiple eQTL with control for population structure and relatedness

PolyQTL is a statistical method to perform multiple eQTL detection, and it has a control for population structure and relatedness with the available genetic relatedness matrix (GRM).

This repository contains source code, and sample data, and the step to run it. If you have any questions or comments, please contact bzeng30@gatech.edu or greg.gibson@biology.gatech.edu

## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation.

### Installation

GCC >=4.7.0, C++ library openmp are needed. 

There is a Makefile, and you can just run "make" to install the package.

#### Running

Genetic relatedness matrix (GRM) is needed to used to control for population structure and relatedness. You can calculate it with external packages, like GCTA, GEMMA.

There are two modes to run the package: 1. Conditional analysis. In this mode, conditional analysis was firstly conducted, and in each iteration, mixed linear model component of GEMMA package was used to detect peak signal. For each detected peak signal, all variants locating in high LD (r2>=0.3) were extracted, and sampling of causal states was run to estimate the importance of explored variants; 2. One-step. In this mode, the step to detect peak signal was skipped, and importance of each variant was estimated.

To run PolyQTL in conditional analysis mode, three files are wanted: 1. GRM; 2. phenotype-genotype file plink bfile format for genotype; 3. binary plink format genotype.

For the phenotype-genotype file, it should be in the format: 

    Ind phe1 phe2 ... phen variant1 variant2 variant3 ... variantM

    Ind1 a11 a12  ... p11 p12 p13 ... p1M

    Ind2 a21 a22  ... p21 p22 p23 ... p2M

    .

    .

    .

    Indn an1 an2   ... pn1 pn2 pn3 ... pnM
 
The package needs some informations to recognise phenotypes, and currently, the phenotype name should be in the format "phe_*". 

To run PolyQTL in one-step mode, three files are needed: 1. GRM; 2. phenotype; 3. genotype dosage.

For one-column phenotype, it should be one-column, and contains the phenotype value

    1.42743
    2.68416
    3.09885
    2.91031
    3.19284
    3.54407
    1.95301
    2.81576
    ......
 
For genotype dosage, it should be in the format below: 
    rs1	rs2	rs3	rs4
    1	2	2	2
    1	1	1	1
    0	2	2	2
    0	2	0	0
    1	2	1	1
    1	2	2	2
    1	2	2	2
    0	1	1	1
    2	2	2	2


##### Examples

You can run the following command to have a sense of how PolyQTL works.

1. Conditional-analysis mode

    ./PolyQTL          -t  1     -P   data/adjust_adjust_combined_CATSPER1_genotype_phenotype  -T  CATSPER1  -G ./GRM_for_CAGE_1763_individual  -o output_test -Z  data/genotype_CATSPER1

2. One-step mode

    ./PolyQTL  -o output_PolyQTL    -p data/CATSPER1.phe  -c 1   -t  1  -x data/CATSPER1.geno -G ./GRM_for_CAGE_1763_individual
