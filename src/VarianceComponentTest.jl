module VarianceComponentTest

# package code goes here
using Distributions
using Docile

include("readgeno!.jl")
include("readAnnotate!.jl")
include("kinshipcoef.jl")
include("vctestnullsim.jl")
include("vctest.jl")
include("loopwinFixsize.jl")
include("loopwinAnnot.jl")
include("gwasvctest.jl")

end # module
