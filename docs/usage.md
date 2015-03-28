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
## Input PLINK files

Option `plinkFile` indicates the file name (without extension) for the input PLINK files. All three PLINK files have the same file name but different extensions. Make sure the three PLINK files are at the same directory.

If the three PLINK files are `plink.bed`, `plink.bim` and `plink.fam`, then use

```julia
gwasvctest(plinkFile = "/PATH/OF/plink")
```

Replace `/PATH/OF/` with the actual path of PLINK files.

---
## Input covariates file

Option `covFile` indicates the file name for the input covariates file. If the covariates file is `covariates.txt`, then use

```julia
gwasvctest(covFile = "/PATH/OF/covariates.txt")
```

If option `covFile` is not specified, the covariates matrix **X** will be automatically set to a *n*-by-1 matrix with all elements equal to 1, where *n* is the number of individuals.

---
## Input trait file

Option `traitFile` indicates the file name for input trait file. If the trait file is `y.txt`, then use

```julia
gwasvctest(traitFile = "/PATH/OF/y.txt")
```

If option `traitFile` is not specified, the response vector **y** will be set automatically to the phenotypes obtained from the **.fam** PLINK file.

---
## Output file

Option `outFile` indicates the file name for the output file. If the output file name is set to `test.out`, then use

```julia
gwasvctest(outFile = "/PATH/OF/test.out")
```

Replace `/PATH/OF/` with the path where you want to store the output file. If option `outFile` is not specified, the output file name will be set to `plinkFile.out` and it will be stored at the current directory.