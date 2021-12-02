module FrictionModels

using DifferentialEquations

export FrictionModel
export CoulombModel
export LuGreModel
export friction

include("FrictionTypes.jl")
include("CoulombModels.jl")

end
