module Adaptive

using WriteVTK

export lglpoints, baryweights, spectralderivative, interpolationmatrix
export creategrid!, computemetric!
export creategrid1d, creategrid2d, creategrid3d, computemetric, creategrid
export writemesh
export nonconinterpolation

include("quadrature.jl")
include("operators.jl")
include("metric.jl")
include("noncon.jl")
include("vtk.jl")


end # module
