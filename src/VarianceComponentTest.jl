module VarianceComponentTest

using Distributions
using Docile

export
       readgeno!,
       readAnnotate!,
       kinshipcoef,
       vctestnullsim,
       vctest,
       loopwinFixsize,
       loopwinAnnot,
       gwasvctest

include("readgeno!.jl")
include("readAnnotate!.jl")
include("kinshipcoef.jl")
include("vctestnullsim.jl")
include("vctest.jl")
include("loopwinFixsize.jl")
include("loopwinAnnot.jl")
include("gwasvctest.jl")

end
