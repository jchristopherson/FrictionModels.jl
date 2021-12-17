# GeneralizedMaxwellSlipModels.jl

function stribeck_curve(
    mdl::GeneralizedMaxwellSlipModel,
    nrm::Number,
    vel::Number
)
    Fc = nrm * mdl.coulomb_coefficient
    Fs = nrm * mdl.static_coefficient
    return Fc + (Fs - Fc) * exp(-(vel / mdl.stribeck_velocity)^2)
end

function element_state_equation(
    mdl::GeneralizedMaxwellSlipModel,
    element::MaxwellElement,
    nrm::Number,
    vel::Number,
    s::Number,
    z::Number
)
    o = one(typeof(nrm))
    si = element.scale_factor * s
    C = element.scale_factor * mdl.attraction_parameter
    if abs(z) <= abs(si)
        dzdt = vel
    else
        dzdt = sign(vel) * C * (o - z / si)
    end
    return dzdt
end

function heuristic_state_equation!(
    mdl::GeneralizedMaxwellSlipModel,
    t::T,
    nrm::T,
    pos::T,
    vel::T,
    z::Array{T},
    dzdt::Array{T}
) where T <: Number

    nelements = min(length(z), length(mdl.elements))
    s = stribeck_curve(mdl, nrm, vel)
    for i in (1:nelements)
        dzdt[i] = element_state_equation(mdl, mdl.elements[i], nrm, vel, s, z[i])
    end
end

function heuristic_force_equation(
    mdl::GeneralizedMaxwellSlipModel,
    t::T,
    nrm::T,
    pos::T,
    vel::T,
    z::Array{T},
    dzdt::Array{T}
) where T <: Number
    
    nelements = min(length(z), length(dzdt), length(mdl.elements))
    F = mdl.viscous_damping * vel
    for i in (1:nelements)
        F = F + mdl.elements[i].stiffness * z[i] +
            mdl.elements[i].damping * dzdt[i]
    end
    return F
end

# ------------------------------------------------------------------------------
function model_from_array(mdl::GeneralizedMaxwellSlipModel, x::Array{T}) where T <: Number
    n = Int(floor(x[1]))
    mu_s = x[2]
    mu_c = x[3]
    c = x[4]
    vs = x[5]
    b = x[6]
    elements = Vector{MaxwellElement}()
    j = 7
    for i in (1:n)
        e = MaxwellElement(x[j], x[j+1], x[j+2])
        j = j + 3
        push!(elements, e)
    end
    return GeneralizedMaxwellSlipModel(elements, mu_s, mu_c, c, vs, b)
end

function model_to_array(x::GeneralizedMaxwellSlipModel)
    c = [
        Float64(length(x.elements)),
        x.static_coefficient,
        x.coulomb_coefficient,
        x.attraction_parameter,
        x.stribeck_velocity,
        x.viscous_damping
    ]
    for i in (1:length(x.elements))
        push!(c, x.elements[i].stiffness)
        push!(c, x.elements[i].damping)
        push!(c, x.elements[i].scale_factor)
    end
    return c
end
