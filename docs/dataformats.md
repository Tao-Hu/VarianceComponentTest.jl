# Data Formats

The type of files involved in this package can be divided into two categories: input files and output file. The summary of different files are listed as

|Name                 |Category |Extension |Description
|:--------------------|:-------:|:--------:|:------------
|Binary genotype file |Input    |.bed      |Binary PLINK file which contains genotypes information
|Marker file          |Input    |.bim      |Ordinary PLINK file which contains information of markers (*e.g.* SNPs)
|Pedigree file        |Input    |.fam      |Ordinary PLINK file which contains information of pedigrees
|Covariates file      |Input    |.txt      |Text file which contains the values of covariates
|Trait file           |Input    |.txt      |Text file which contains the values of traits/phenotypes
|Kinship file         |Input    |.txt      |Text file which contains the values of kinship matrix
|Output file          |Output   |.out      |Flat file which contains p-values for each group of markers under certain testing scheme

---
## PLINK files

Three types of [PLINK files](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed) are needed: **.bed**, **.bim** and **.fam** files.

### BED file

BED file is a binary PLINK file which contains genotypes information. More information about the format of BED file is [here](http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml).

**Note**: if you try to open the BED file in a normal editor, you will see unreadable characters.

The data stored in BED file are in units of byte. The actual genotypes are beginning from the fourth byte since the first 3 bytes have special meanings. Here, what we concern is the third byte which indicates whether the BED file is in SNP-major or individual-major mode: a value of 00000001 indicates SNP-major, and a value of 00000000 indicates individual-major. Only the SNP-major mode is accepted in our program, otherwise an error will be given.

### BIM file

BIM file is an ordinary PLINK file which contains information of markers (*e.g.* SNPs). In BIM file, each line represents a SNP. The fields in BIM file are

* Chromosome
* SNP ID
* Genetic distance (the unit is Morgan)
* Physical position (the unit is bp)
* Allele 1
* Allele 2

Fields are separated by a whitespace. A typical BIM file looks like

```
3	3-47892183	0	47892183	T	C
3	3-47892350	0	47892350	C	T
3	3-47892383	0	47892383	A	G
3	3-47892431	0	47892431	T	G
3	3-47892574	0	47892574	T	A
...
```

### FAM file

FAM file is an ordinary PLINK file which contains information of pedigrees. In FAM file, each line represents an individual. The fields in FAM file are

* Family ID
* Individual ID
* Paternal ID (0 means that individual is a founder)
* Maternal ID (0 means that individual is a founder)
* Sex (1 = male; 2 = female)
* Phenotype

Fields are separated by a whitespace. A typical FAM file looks like

```
2 T2DG0200018 0 0 1 -9
2 T2DG0200023 0 0 2 -9
2 T2DG0200027 0 0 2 -9
2 T2DG0200031 T2DG0200001 T2DG0200015 1 -9
2 T2DG0200032 T2DG0200001 T2DG0200015 2 -9
...
```

---
## Covariates file

Covariates file is a text file which contains the values of covariates. In covariates file, each line represents an individual. The first two fields are `Family ID` and `Individual ID`, the rest fields are the covariates/predictors (*e.g.* `sex`, `age`, `weight` and *etc*). The intercept should NOT be included. Fields are separated by a whitespace. A typical covariates file look like

```
2.0 T2DG0200001 1.0 0.8333876751238651 0.4668769615832563 -0.7067078743064089 0.4914401659107072
2.0 T2DG0200002 1.0 1.560273019093175 1.321338604534663 -0.31604849390507345 -0.25920628775255705
2.0 T2DG0200003 1.0 -1.421681558701967 0.5815057104664987 -0.8991654720299747 0.5865043498666274
2.0 T2DG0200004 0.0 -0.7354246339009091 1.3591699397691537 0.6110166359709753 0.5707197234022439
2.0 T2DG0200005 0.0 0.4386793918374321 1.2766798171425295 -0.7428555065380336 -0.03767993395934451
...
```

If no covariates file is provided, the covariates matrix **X** will be automatically set to a *n*-by-1 matrix with all elements equal to 1, where *n* is the number of individuals.

---
## Trait file

Trait file is a text file which contains the values of traits/phenotypes. In trait file, each line represents an individual. The fields in trait file are

* Family ID
* Individual ID
* Trait value

Fields are separated by a whitespace. A typical trait file looks like

```
2.0 T2DG0200001 5.0860722020317315
2.0 T2DG0200002 3.621930716051974
2.0 T2DG0200003 0.7356384799994807
2.0 T2DG0200004 2.8070180572616326
2.0 T2DG0200005 2.735167495689626
...
```

If no trait file is provided, the response vector **y** will be set automatically to the `Phenotype` field obtained from the FAM file.

---
## Kinship file
Kinship file is a text file which contains the values of kinship matrix. If the data set which you are going to analysis contains *n* individuals, then the dimension of kinship matrix is *n*-by-*n*, which means there should be *n* rows and *n* fields in the kinship file, and fields are separated by comma.

For example, if a data set contains 5 individuals, the kinship file looks like

```
0.5,0,0,0.25,0
0,0.5,0,0,0
0,0,0.5,0,0
0.25,0,0,0.5,0
0,0,0,0,0.5
```

---
## Output file

Output file is a flat file which contains p-values for each group of markers under certain testing scheme. In output file, each line represents the p-value for one group of markers (the first line is the header). The fields in output file are

* Starting marker ID of group
* Ending marker ID of group
* p-value for group

Fields are separated by comma. A typical output file looks like

```
StartSNP,EndSNP,pvalue
3-47892183,3-47905079,1.0
3-47905164,3-47920240,1.0
3-47920887,3-47931984,1.0
3-47932222,3-47946708,1.0
3-47946801,3-47959188,0.08519730073277564
...
```