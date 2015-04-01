# Usage

The core function in *ExactVarianceComponentTest.jl* is *gwasvctest*, which wraps the three types of exact tests. The usage for *gwasvctest* is

```julia
gwasvctest(Name = Value)
```

Where `Name` is the option name and `Value` is the value of the corresponding option. Note that it's okay if `Value` is not provided since every option has a default value.

There are two ways to call *gwasvctest*

* Open up Julia and type

```julia
julia> using ExactVarianceComponentTest
julia> gwasvctest(Name = Value)
```

* From command line, type

```
$ julia -E 'using ExactVarianceComponentTest; gwasvctest(Name = Value)'
```

Next, all options are listed below.

---
## Background

??? Need to add some background introductions for the algorithm.

---
## Specify input PLINK files

Option `plinkFile` indicates the file name (without extension) for the input PLINK files. All three PLINK files have the same file name but different extensions. Make sure the three PLINK files are at the same directory.

If the three PLINK files are *plink.bed*, *plink.bim* and *plink.fam*, then use

```julia
gwasvctest(plinkFile = "/PATH/OF/plink")
```

Replace "/PATH/OF/" with the actual path of PLINK files.

---
## Specify input covariates file

Option `covFile` indicates the file name for the input covariates file. If the covariates file is *covariates.txt*, then use

```julia
gwasvctest(covFile = "/PATH/OF/covariates.txt")
```

If option `covFile` is not specified, the covariates matrix **X** will be automatically set to a *n*-by-1 matrix with all elements equal to 1, where *n* is the number of individuals.

---
## Specify input trait file

Option `traitFile` indicates the file name for input trait file. If the trait file is *y.txt*, then use

```julia
gwasvctest(traitFile = "/PATH/OF/y.txt")
```

If option `traitFile` is not specified, the response vector **y** will be set automatically to the phenotypes obtained from the **.fam** PLINK file.

---
## Specify output file

Option `outFile` indicates the file name for the output file. If the output file name is set to *test.out*, then use

```julia
gwasvctest(outFile = "/PATH/OF/test.out")
```

Replace "/PATH/OF/" with the path where you want to store the output file. If option `outFile` is not specified, the output file name will be set to *plinkFile-julia.out* and it will be stored at the current directory.

---
## Choose testing scheme

*ExactVarianceComponentTest.jl* provides three types of exact tests: exact likelihood ratio test (eLRT), exact restricted likelihood ratio test (eRLRT) and exact score test (eSC). Option `test` indicates which testing scheme you want to perform. The usage is

* `gwasvctest(test = "eLRT")`: perform exact likelihood ratio test
* `gwasvctest(test = "eRLRT")`: perform exact restricted likelihood ratio test
* `gwasvctest(test = "eScore")`: perform exact score test

The default value for option `test` is *eRLRT*.

---
## Choose method for obtaining null distribution of test statistic and computing p-value

Option `pvalueComputing` indicates which method will be used to obtain the null distribution of test statistic. The usage is

* `gwasvctest(pvalueComputing = "MonteCarlo")`: use a Monte Carlo method by generating many replicates to obtain the exact null distribution of test statistic
* `gwasvctest(pvalueComputing = "chi2")`: use a mixed Chi squared distribution to approximate the null distribution of test statistic. Only valid under **eLRT** or **eRLRT**

The default value for option `pvalueComputing` is *chi2*. The approximation effect of mixed Chi squared distribution has been showed to be good enough, and such approximation will be faster than Monte Carlo method since much less replicates need to be generated. So please use *chi2* whenever possible.

---
## Choose number of replicates to generate for obtaining null distribution of test statistic

Option `nNullSimPts` lets you to decide how many replicates to generate for obtaining null distribution of test statistic.

* For `pvalueComputing = "MonteCarlo"`, the more replicates to generate, the more precise the p-value will be (smaller standard error)
* For `pvalueComputing = "chi2"`, the number of replicates does not count too much, it only effect the estimate of probability at zero for test statistic

`nNullSimPts` should take positive integer, and the default value is 10,000.

---
## Choose method for computing kinship matrix

Option `kinship` indicates how to obtain the kinship matrix. The usage is

* `gwasvctest(kinship = "GRM")`: compute kinship matrix by genetic relationship matrix (GRM)
* `gwasvctest(kinship = "MoM")`: compute kinship matrix by method of moments (MoM)
* `gwasvctest(kinship = "theoretical")`: compute kinship matrix theoretically
* `gwasvctest(kinship = "none")`: replace *none* with the file name in which store your own pre-calculated kinship matrix

The default value for option `kinship` is *GRM*.

---
## Choose method for computing kernel matrix

Option `kernel` indicates how to obtain kernel matrix. The usage is

* `gwasvctest(kernel = "GRM")`: compute kernel matrix by genetic relationship matrix (GRM)
* `gwasvctest(kernel = "IBS1")`
* `gwasvctest(kernel = "IBS2")`
* `gwasvctest(kernel = "IBS3")`

The default value for option `kernel` is *GRM*.

---
## Determine the rank used to approximate the kinship matrix

In our algorithm, we do a low rank approximation of kinship matrix to reduce the multiple variances components testing problem to two variance components problem. However, there is a trade-off for choosing the rank for approximating: if the rank is high, then a better approximation can be obtained, but a high rank will lead to a small signal-to-noise ratio, which will decrease the power of test.

Option `infLambda` provides a way to control such trade-off. Theoretically, `infLambda` can take any real value. If `infLambda` takes non positive value, then the algorithm will take the highest possible approximating rank. If `infLambda` takes a positive value, the algorithm will take a smaller approximating rank, and the larger the value is, the smaller the approximating rank will be. The default value is 0.

For example, want to take the highest rank, then

```julia
gwasvctest(infLambda = 0.0)
```

---
## Determine the group size for each SNP set

Option `windowSize` determines the group size of a SNP set for which you want to test. It can only take positive integer. The default value is 50.

For example, want to set the group size to 60, then

```julia
gwasvctest(windowSize = 60)
```