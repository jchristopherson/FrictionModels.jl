module FrictionModels

using LsqFit
using Interpolations
using DifferentialEquations

export FrictionModel
export CoulombModel
export LuGreModel
export HyperbolicModel
export friction
export fit_model

include("FrictionTypes.jl")
include("HelperRoutines.jl")
include("Fitting.jl")
include("CoulombModels.jl")
include("LuGreModels.jl")
include("HyperbolicModels.jl")

end
