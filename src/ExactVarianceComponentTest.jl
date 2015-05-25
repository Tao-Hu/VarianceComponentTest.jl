module ExactVarianceComponentTest

# package code goes here
using Distributions

include("readgeno!.jl")
include("readAnnotate!.jl")
include("kinshipcoef.jl")
include("vctestnullsim.jl")
include("vctest.jl")
include("loopwinFixsize.jl")
include("loopwinAnnot.jl")
include("gwasvctest.jl")

end # module
