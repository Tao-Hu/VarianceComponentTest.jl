# Data Formats

The type of files involved in this package can be divided into two categories: input files and output file. The summary of different files are listed as

|Name                 |Category |Extension |Description
|:--------------------|:-------:|:--------:|:------------
|Binary genotype file |Input    |.bed      |Binary PLINK file that contains genotypes information
|Marker file          |Input    |.bim      |Ordinary PLINK file that contains information of markers (*e.g.* SNPs)
|Pedigree file        |Input    |.fam      |Ordinary PLINK file that contains information of pedigrees
|Covariates file      |Input    |.txt      |Text file that contains the values of covariates
|Trait file           |Input    |.txt      |Text file that contains the values of traits/phenotypes
|Annotation file      |Input    |.txt      |Text file that contains the annotation information of markers
|Kinship file         |Input    |.txt      |Text file that contains the values of kinship matrix
|Output file          |Output   |.out      |Flat file that contains p-values for each group of markers under certain testing scheme

---
## PLINK files

Three types of [PLINK files](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed) are needed: **.bed**, **.bim** and **.fam** files.

### BED file

BED file is a binary PLINK file that contains genotypes information. More information about the format of BED file is [here](http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml).

**Note**: if you try to open the BED file in a normal editor, you will see unreadable characters.

The data stored in BED file are in units of byte. The actual genotypes are beginning from the fourth byte since the first 3 bytes have special meanings. Here, what we concern is the third byte which indicates whether the BED file is in SNP-major or individual-major mode: a value of 00000001 indicates SNP-major, and a value of 00000000 indicates individual-major. Only the SNP-major mode is accepted in our program, otherwise an error will be given.

### BIM file

BIM file is an ordinary PLINK file that contains information of markers (*e.g.* SNPs). In BIM file, each line represents a SNP. The fields in BIM file are

* Chromosome code (either an integer, or "X"/"Y"/"XY"/"MT", "0" indicates unknown)
* SNP ID
* Genetic distance (the unit is Morgan)
* Physical position (the unit is bp)
* Allele 1
* Allele 2

Make sure the markers are in increasing order in terms of `Chromosome` and `Physical position` (This is usually the case for a BIM file).

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

FAM file is an ordinary PLINK file that contains information of pedigrees. In FAM file, each line represents an individual. The fields in FAM file are

* Family ID
* Individual ID
* Paternal ID (0 means that individual is a founder)
* Maternal ID (0 means that individual is a founder)
* Sex (1 = male; 2 = female; 0 = unknown)
* Phenotype

For the case-control data, the `Phenotype` field is defined as: "1" = control, "2" = case, "-9"/"0"/non-numeric = missing data.

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

Covariates file is a text file that contains the values of covariates. In covariates file, each line represents an individual. The format of covariates file is the same as that in PLINK: the first two fields are `Family ID` and `Individual ID`, the rest fields are the covariates/predictors (*e.g.* *sex*, *age*, *weight* and *etc*). The intercept should NOT be included. Fields are separated by a whitespace. A typical covariates file look like

```
2 T2DG0200001 1.0 0.8333876751238651 0.4668769615832563 -0.7067078743064089 0.4914401659107072
2 T2DG0200002 1.0 1.560273019093175 1.321338604534663 -0.31604849390507345 -0.25920628775255705
2 T2DG0200003 1.0 -1.421681558701967 0.5815057104664987 -0.8991654720299747 0.5865043498666274
2 T2DG0200004 0.0 -0.7354246339009091 1.3591699397691537 0.6110166359709753 0.5707197234022439
2 T2DG0200005 0.0 0.4386793918374321 1.2766798171425295 -0.7428555065380336 -0.03767993395934451
...
```

If a specific covariate value is missing, please use "NaN" to indicate that missing value

If no covariates file is provided, the covariates matrix **X** will be automatically set to a *n*-by-1 matrix with all elements equal to 1, where *n* is the number of individuals.

---
## Trait file

Trait file is a text file that contains the values of traits/phenotypes. In trait file, each line represents an individual. The format of trait file is the same as phenotype file in PLINK. The fields in trait file are

* Family ID
* Individual ID
* Trait value ("NaN" for missing value)

Fields are separated by a whitespace. A typical trait file looks like

```
2 T2DG0200001 5.0860722020317315
2 T2DG0200002 3.621930716051974
2 T2DG0200003 0.7356384799994807
2 T2DG0200004 2.8070180572616326
2 T2DG0200005 2.735167495689626
...
```

If no trait file is provided, the response vector **y** will be set automatically to the `Phenotype` field obtained from the FAM file.

---
## Annotation file
Annotation file is a text file that contains the information about in which gene a SNP belongs to. The annotation file must satisfy two conditions: (1) it should contain annotation information of all markers in the BIM file, (2) the annotated markers in the annotation file should in increasing order in terms of which chromosome a marker is located and its physical position.

In annotation file, each line represents a marker and its annotation information. The fields in annotation file are

* Gene name (indicate which gene the SNP belongs to)
* SNP ID

Fields are separated by a comma. A typical trait file looks like

```
gene1,3-47892183
gene1,3-47892350
gene1,3-47892383
gene2,3-47892431
gene2,3-47892574
...
```

---
## Kinship file
Kinship file is a text file that contains the values of kinship matrix. If the data set which you are going to analysis contains *n* individuals, then the dimension of kinship matrix is *n*-by-*n*, which means there should be *n* rows and *n* fields in the kinship file, and fields are separated by comma.

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

Output file is a flat file that contains p-values and other information for each group of markers under certain testing scheme. In output file, each line represents for one group of markers (the first line is the header). There are two formats of the output file: one is for the case when no annotation file is provided, the other is for the case when an annotation file is provided.

#### No annotation file provided

The fields in output file are

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

#### Annotation file provided

The fields in output file are

* Starting marker ID of gene
* Ending marker ID of gene
* Chromosome of gene
* Starting physical position of markers in the gene
* Ending physical position of markers in the gene
* Gene name
* Number of markers contains in the gene
* p-value for gene

Fields are separated by comma. A typical output file looks like

```
StartSNP,EndSNP,Chr,StartPos,EndPos,GeneName,#SNPs,pvalue
1-865584,1-879184,1,865584,879184,SAMD11,8,0.36926357377277497
1-880922,1-892388,1,880922,892388,NOC2L,9,0.4102919217592677
1-897009,1-900371,1,897009,900371,KLHL17,10,0.2068731433962015
rs62639980,1-909723,1,901922,909723,PLEKHN1,13,0.19108809387378453
1-934982,1-934982,1,934982,934982,HES4,1,0.31590801662049744
...
```