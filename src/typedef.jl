# This file define the structures and subfunctions used in the program

# type of stats
type Stats
  vc1_pvalue::Float64
  vc1_pvalue_se::Float64
  vc1_teststat::Float64
  iterations::Int64
  method::String
  residual::Array{Float64,1}
  test::String
end

Stats() = Stats(0.0, 0.0, 0.0, 10, "MLE", Float64[], "eLRT")

# type of V in the form of eigen decomposition
type eigenV
  U
  eval
end

eigenV() = eigenV([], [])
