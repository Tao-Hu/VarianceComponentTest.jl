# Variance Component Test

*VarianceComponentTest.jl* is a Julia package for performing exact variance component tests in genome-wide association study (GWAS). It provides three types of exact tests

* exact likelihood ratio test (*eLRT*)
* exact restricted likelihood ratio test (*eRLRT*)
* exact score test (*eScore*)

The input files for *VarianceComponentTest.jl* are PLINK formatted files (**.bed**, **.bim** and **.fam** file), covariates file (**.txt** file) and trait file (**.txt** file). You should have these input files prepared before running our program. The output file (**.out** file) is a simple comma-delimited file containing the p-values for each group of SNPs under a certain testing scheme (eLRT, eRLRT or eScore).

To use *VarianceComponentTest.jl*, you need to call the *gwasvctest()* function.

## Contents

* [Installation](installation.md)
 
* [Dataformats](dataformats.md)

* [Usage](usage.md)

* [Examples](./examples/examples.md)
