module FrictionModels

using LsqFit
using Interpolations
using DifferentialEquations

export FrictionModel
export HeuristicFrictionModel
export CoulombModel
export LuGreModel
export HyperbolicModel
export MaxwellElement
export GeneralizedMaxwellSlipModel
export friction
export fit_model

include("FrictionTypes.jl")
include("HelperRoutines.jl")
include("HeuristicFrictionModels.jl")
include("Fitting.jl")
include("CoulombModels.jl")
include("LuGreModels.jl")
include("HyperbolicModels.jl")
include("GeneralizedMaxwellSlipModels.jl")

Base.length(::FrictionModel) = 1
Base.size(::FrictionModel) = ()
Base.getindex(x::FrictionModel, i) = x
Base.iterate(x::FrictionModel, state = 1) = state > length(x) ? nothing : (x, state = state + 1)

end
