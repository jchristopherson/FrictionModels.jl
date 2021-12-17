# FrictionTypes.jl

"""
Abstract type defining a basic friction model.
"""
abstract type FrictionModel end

"""
An abstract heuristic friction model.
"""
abstract type HeuristicFrictionModel <: FrictionModel end

"""
A simple Coulomb friction model the represents friction force proportional
to the normal force via a single friction coefficient.
"""
struct CoulombModel <: FrictionModel
    """
    The Coulomb friction coefficient.
    """
    coefficient::Number
end

"""
The Lu-Gre friction model.
"""
struct LuGreModel <: HeuristicFrictionModel
    static_coefficient::Number
    coulomb_coefficient::Number
    stribeck_velocity::Number
    bristle_stiffness::Number
    bristle_damping::Number
    viscous_damping::Number
end

"""
A friction model based upon hyperbolic functions proposed by Rodriguez et. al.
"""
struct HyperbolicModel <: FrictionModel
    friction_coefficient::Number
    normalization_coefficient::Number
    dissipation_coefficient::Number
    hysteresis_coefficient::Number
    stribeck_velocity::Number
    viscous_damping::Number
end

"""
A single Maxwell type element used in Maxwell-Slip models.
"""
struct MaxwellElement
    stiffness::Number
    damping::Number
    scale_factor::Number
end

"""
The Generalized Maxwell-Slip Model
"""
struct GeneralizedMaxwellSlipModel <: HeuristicFrictionModel
    elements::Array{MaxwellElement}
    static_coefficient::Number
    coulomb_coefficient::Number
    attraction_parameter::Number
    stribeck_velocity::Number
    viscous_damping::Number
end
