module FrictionModels

using LsqFit
using Interpolations
using DifferentialEquations

export FrictionModel
export CoulombModel
export LuGreModel
export friction
export fit_model

include("FrictionTypes.jl")
include("CoulombModels.jl")
include("LuGreModels.jl")

end
