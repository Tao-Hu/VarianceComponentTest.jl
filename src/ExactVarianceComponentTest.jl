module ExactVarianceComponentTest

# package code goes here
using Distributions

include("typedef.jl")
include("readgeno!.jl")
include("readAnnotate!.jl")
include("kinshipcoef.jl")
include("vctestnullsim.jl")
include("vctest.jl")
include("gwasvctest.jl")

end # module
