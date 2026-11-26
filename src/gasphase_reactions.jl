"""
    get_twobody_rate_constant(α, β, γ, T_gas)

Get the rate constant of a gas-phase twobody reaction,
using the Arrhenius-Kooij equation.

``k = \\alpha (T_{\\mathrm{gas}}/300)^\\beta e^{-\\gamma/T_{\\mathrm{gas}}}``"""
function get_twobody_rate_constant(α, β, γ, T_gas)
    return α * (T_gas / 300.0)^β * get_arrhenius_probability(γ, T_gas)
end

"""Equation 1 from Wakelam et al, 2024."""
function get_ionopol1_rate_constant(α, β, γ, T_gas)
    return α * β * (0.62 + 0.4767 * γ * √(300.0 / T_gas))
end

"""Equation 2 from Wakelam et al, 2024."""
function get_ionopol2_rate_constant(α, β, γ, T_gas)
    return α * β * (1.0 + 0.0967 * γ * √(300.0 / T_gas) + (γ^2 / 10.526) * (300.0 / T_gas))
end

"""
    get_photo_rate_constant(α, γ, radfield, A_V)

Get the rate constant of a photo-reaction.

``k = \\alpha * F_{\\mathrm{UV}} e^{-\\gamma * A_{\\mathrm{V}}} / 1.7``"""
function get_photo_rate_constant(α, γ, radfield, A_V)
    return α * radfield * exp(-γ * A_V) / 1.7
end

"""
    get_crp_photo_rate_constant(α, β, γ, ω_dust, ζ, T_gas)

Get the rate constant of a CRPHOTO reaction.

``k = α * γ / (1-ω_{\\mathrm{dust}}) * ζ * (T_{\\mathrm{gas}}/300)^{\\beta}``
"""
function get_crp_photo_rate_constant(α, β, γ, ω_dust, ζ, T_gas)
    return α * γ / (1.0 - ω_dust) * ζ * (T_gas / 300.0)^β
end

"""
    get_crp_rate_constant(α, ζ)

``k = αζ``
"""
function get_crp_rate_constant(α, ζ)
    return α * ζ
end

@register_symbolic get_twobody_rate_constant(α, β, γ, T_gas)
@register_symbolic get_ionopol1_rate_constant(α, β, γ, T_gas)
@register_symbolic get_ionopol2_rate_constant(α, β, γ, T_gas)
@register_symbolic get_crp_photo_rate_constant(α, β, γ, ω_dust, ζ, T_gas)
@register_symbolic get_crp_rate_constant(α, ζ)
