# Examples

This example illustrates the usages of different options by analyzing a small data set obtained from [GAW18](http://www.gaworkshop.org/gaw18/index.html). This demo data set contains 894 SNPs and 849 individuals. The required/optional input files are

* PLINK file: `chr3-geno-MAP4-849.bed`, `chr3-geno-MAP4-849.bim` and `chr3-geno-MAP4-849.fam`
* Covariates file: `covariates-julia.txt` which contains 5 covariates and intercept
* Trait file: `y-julia.txt`
* Annotation file: `annotation.txt`
* Kinship file: `gaw18_849_kinship.txt`

These data files come with our package, and they are available at [*here*](https://github.com/Tao-Hu/VarianceComponentTest.jl/tree/master/docs/examples).

---
## Basic usage

Three types of exact tests can be performed. Open up a Julia session and type

* exact likelihood ratio test (eLRT)

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eLRT")
```

Then the output will be written to `chr3-geno-MAP4-849-julia.out` at the same directory as PLINK files.

* exact restricted likelihood ratio test (eRLRT)

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eRLRT")
```

* exact score test (eScore)

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eScore")
```

**Note**: if the input files are not at the current directory, you should specify the paths correctly.

You can also call *gwasvctest* from command line. For example, to perform eRLRT

```
$ julia -E 'using VarianceComponentTest; gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eRLRT")'
```

---
## Option `pvalueComputing`

Chi squared approximation is recommended (though you don't have to write it out specifically)

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eRLRT", pvalueComputing = "chi2")
```

If you want to use Monte Carlo method

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eRLRT", pvalueComputing = "MonteCarlo")
```

**Note**: Option `pvalueComputing` is only for eLRT and eRLRT. For eScore, it employs the method of inverting characteristics function to compute the p-value.

---
## Option `nNullSimPts`

If Monte Carlo method is chosen, generate 10,000 replicates should be sufficient (you can also generate even more replicates to increase the precision of p-value)

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eRLRT", pvalueComputing = "MonteCarlo", nNullSimPts = 10000)
```

If Chi squared approximation method is chosen, you can generate more replicates (default value is 1,000,000) to have a better estimator of point mass at 0

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eRLRT", pvalueComputing = "chi2", nNullSimPts = 10000000)
```

---
## Option `annotationFile`

If the annotation information of the markers is available, you can share it by inputting an annotation file with option `annotationFile`

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eRLRT", annotationFile = "annotation.txt")
```

---
## Option `kinship`

For the analysis of unrelated data, the kinship matrix should not be included in the model. So you should specify option `kinship` as *none*

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", test = "eRLRT", kinship = "none")
```

If you have family data and have a file (*gaw18_849_kinship.txt* in this example) which contains the kinship matrix

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", test = "eRLRT", kinship = "gaw18_849_kinship.txt")
```

It's okay if you don't have the kinship file

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", test = "eRLRT", kinship = "GRM")
```

---
## Option `infLambda`

If you are concerning about too high approximating rank for the kinship matrix, use option `infLambda` to set a lower rank

```julia
julia> using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eRLRT", infLambda = 1.0)
```

**Note**: If you are analyzing unrelated data, you can ignore option `infLambda` since the kinship matrix will not be included in the model in unrelated case.

---
## Parallel computing

If your machine has multiple cores, you can active the parallel computing mode to speed up the calculations. For example, you want to distribute the computings to 4 processes

```julia
julia> addprocs(4)
julia> @everywhere using VarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eRLRT", annotationFile = "annotation.txt")
```

You can also call the function from command line

```
$ julia -p 4 -E '@everywhere using VarianceComponentTest; gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "none", test = "eRLRT", annotationFile = "annotation.txt")'
```
