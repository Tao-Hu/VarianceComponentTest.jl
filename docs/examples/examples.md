# Examples

This example illustrates the usages of different options by analyzing a small data set obtained from [GAW18](http://www.gaworkshop.org/gaw18/index.html). This demo data set contains 894 SNPs and 849 individuals. The required/optional input files are

* PLINK file: `chr3-geno-MAP4-849.bed`, `chr3-geno-MAP4-849.bim` and `chr3-geno-MAP4-849.fam`
* Covariates file: `covariates-julia.txt` which contains 5 covariates and intercept
* Trait file: `y-julia.txt`
* Kinship file: `gaw18_849_kinship.txt`

These data files come with our package.

---
## Basic usage

Three types of exact tests can be performed. Open up a Julia session and type

* exact likelihood ratio test (eLRT)

```julia
julia> using ExactVarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eLRT")
```

Then the output will be written to `chr3-geno-MAP4-849.out` at current directory.

* exact restricted likelihood ratio test (eRLRT)

```julia
julia> using ExactVarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eRLRT")
```

* exact score test (eSC)

```julia
julia> using ExactVarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eScore", pvalueComputing = "MonteCarlo")
```

**Note**: 1) option `pvalueComputing` must be *MonteCarlo* under eSC. 2) if the input files are not at the current directory, you should specify the paths correctly.

You can also call *gwasvctest* from command line. For example, to perform eRLRT

```
$ julia -E 'using ExactVarianceComponentTest; gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eRLRT")'
```

---
## Option `pvalueComputing`

For eLRT and eRLRT, Chi squared approximation is recommended (though you don't have to write it out specifically)

```julia
julia> using ExactVarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eRLRT", pvalueComputing = "chi2")
```

For eSC, Monte Carlo method is the only option

```julia
julia> using ExactVarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eScore", pvalueComputing = "MonteCarlo")
```

---
## Option `nNullSimPts`

If Monte Carlo method is chosen and you want to increase the precision of p-value, then you can generate more replicates

```julia
julia> using ExactVarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eScore", pvalueComputing = "MonteCarlo", nNullSimPts = 100000)
```

---
## Option `kinship`

It's okay if you don't have the kinship file

```julia
julia> using ExactVarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eRLRT", kinship = "GRM")
```

---
## Option `infLambda`

If you are concerning about too high approximating rank, use option `infLambda` to set a lower rank

```julia
julia> using ExactVarianceComponentTest
julia> gwasvctest(plinkFile = "chr3-geno-MAP4-849", covFile = "covariates-julia.txt", traitFile = "y-julia.txt", kinship = "gaw18_849_kinship.txt", test = "eRLRT", infLambda = 1.0)
```