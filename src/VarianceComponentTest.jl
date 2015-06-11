module VarianceComponentTest

import Distributions
if VERSION < v"0.4-"
       using Docile
end

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
