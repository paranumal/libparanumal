module Adaptive

export lglpoints, baryweights, spectralderivative, interpolationmatrix
export creategrid!, computemetric!
export creategrid1d, creategrid2d, creategrid3d, computemetric, creategrid

include("quadrature.jl")
include("operators.jl")
include("metric.jl")


end # module
