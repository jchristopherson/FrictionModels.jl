# FrictionTypes.jl

"""
Abstract type defining a basic friction model.
"""
abstract type FrictionModel end

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

struct LuGreModel <: FrictionModel
    static_coefficient::Number
    coulomb_coefficient::Number
    stribeck_velocity::Number
    bristle_stiffness::Number
    bristle_damping::Number
    viscous_damping::Number
end