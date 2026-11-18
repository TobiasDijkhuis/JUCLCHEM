const H2_FLUX_CT = 1.0e-10

"""
    get_thermal_velocity(mass, T)


``v_{\\mathrm{th}}=\\sqrt{\\frac{8k_\\mathrm{B}T}{πm}}``
"""
function get_thermal_velocity(mass::Real, T::Real)::Real
    thermal_vel = √(8 * K_B_SI * T / (π * AMU_SI * mass))
    return thermal_vel
end

"""
    get_arrhenius_probability(E_reac, temperature)

``P = e^{-E/T}``"""
function get_arrhenius_probability(E_reac::Real, temperature::Real)::Real
    return exp(-E_reac / temperature)
end

"""
    get_tunneling_probability(E_reac, μ_ij, barrier_width)

``P = e^{\\frac{-2a}{\\hbar}\\sqrt{2\\mu_{ij}E_{\\mathrm{reac}}}}``"""
function get_tunneling_probability(E_reac::Real, μ_ij::Real, barrier_width::Real)::Real
    return exp(-2 * barrier_width * 1.0e-10 / ħ_SI * √(2 * μ_ij * AMU_SI * E_reac * K_B_SI))
end

"""
    get_reaction_probability(E_reac, μ_ij, T_dust, barrier_width=1.4)

Get the total reaction probability (including tunneling).

Maximum of [`get_arrhenius_probability`](@ref) and [`get_tunneling_probability`](@ref)
"""
function get_reaction_probability(
    E_reac::Real,
    μ_ij::Real,
    T_dust::Real,
    barrier_width::Real = 1.4,
)::Real
    return max(
        get_tunneling_probability(E_reac, μ_ij, barrier_width),
        get_arrhenius_probability(E_reac, T_dust),
    )
end

"""
    get_sticking_coefficient_Hollenbach_McKee(T_gas, T_dust)

Equation 3.7 from Hollenbach & McKee 1979.
"""
function get_sticking_coefficient_Hollenbach_McKee(T_gas::Real, T_dust::Real)::Real
    S = 1 / (1 + 0.04 * √(T_gas + T_dust) + 0.2 * T_gas / 100 + 0.08 * (T_gas / 100)^2)
    return S
end

"""
    _get_factor_1(μ, ν, E_H2, T_dust)

The first term in equation 15 from Cazaux & Tielens, 2002.
"""
function _get_factor_1(μ::Real, ν::Real, E_H2::Real, T_dust::Real)::Real
    return μ * H2_FLUX_CT / (2 * ν * get_arrhenius_probability(E_H2, T_dust))
end

"""
    _get_factor_2(E_HC, E_S, E_HP, T_dust)

Equation 17 from Cazaux & Tielens, 2002."""
function _get_factor_2(E_HC::Real, E_S::Real, E_HP::Real, T_dust::Real)::Real
    return (1 + √((E_HC - E_S) / (E_HP - E_S)))^2 /
        4 *
        get_arrhenius_probability(E_S, T_dust)
end

"""
    _get_ϵ(ν, E_HC, E_S, E_HP, T_dust)

Equation 16 from Cazaux & Tielens, 2002.
"""
function _get_ϵ(ν::Real, E_HC::Real, E_S::Real, E_HP::Real, T_dust::Real)::Real
    return 1 /
        (
            1 +
                ν /
                    (2 * H2_FLUX_CT) *
                    get_arrhenius_probability(-1.5 * E_HC, T_dust) *
                    (1 + √((E_HC - E_S) / (E_HP - E_S)))^2
        )
end

const GRAIN_DATA = Base.ImmutableDict(
    "Si" => Base.ImmutableDict(
        "μ" => 0.005,
        "E_S" => 110.0,
        "E_H2" => 320.0,
        "E_HP" => 450.0,
        "E_HC" => 3.0e4,
        "ν_H2" => 3.0e12,
        "ν_H" => 1.3e13,
        "CROSS_SECTION" => 8.437e-26,
    ),
    "Gr" => Base.ImmutableDict(
        "μ" => 0.005,
        "E_S" => 260.0,
        "E_H2" => 520.0,
        "E_HP" => 800.0,
        "E_HC" => 3.0e4,
        "ν_H2" => 3.0e12,
        "ν_H" => 1.3e13,
        "CROSS_SECTION" => 7.908e-26,
    ),
)

function get_formation_efficiency(material::AbstractString, T_dust::Real)::Real
    μ = GRAIN_DATA[material]["μ"]
    E_H2 = GRAIN_DATA[material]["E_H2"]
    E_HC = GRAIN_DATA[material]["E_HC"]
    E_S = GRAIN_DATA[material]["E_S"]
    E_HP = GRAIN_DATA[material]["E_HP"]
    ν_H2 = GRAIN_DATA[material]["ν_H2"]
    ν_H = GRAIN_DATA[material]["ν_H"]
    factor_1 = _get_factor_1(μ, ν_H2, E_H2, T_dust)
    factor_2 = _get_factor_2(E_HC, E_S, E_HP, T_dust)
    ϵ = _get_ϵ(ν_H, E_HC, E_S, E_HP, T_dust)
    return ϵ / (1 + factor_1 + factor_2)
end

function get_CT_H2_form_efficiency(T_gas::Real, T_dust::Real)::Real
    v_th = get_thermal_velocity(AtomicMassConstant, T_gas)
    sticking_coefficient = get_sticking_coefficient_Hollenbach_McKee(T_gas, T_dust)

    graphite_efficiency = get_formation_efficiency("Gr", T_dust)
    graphite_rate = 0.5 *
        v_th *
        sticking_coefficient *
        GRAIN_DATA["Gr"]["CROSS_SECTION"] *
        graphite_efficiency
    silicon_efficiency = get_formation_efficiency("Si", T_dust)
    silicon_rate = 0.5 *
        v_th *
        sticking_coefficient *
        GRAIN_DATA["Si"]["CROSS_SECTION"] *
        silicon_efficiency

    return graphite_rate + silicon_rate
end

"""
    get_ν_HH(E_bind, mass)

``\\nu_{\\mathrm{HH}}=\\sqrt{\\frac{2*SURFACE_SITE_DENSITY*E_{\\mathrm{bind}}}{π^2m}}``
"""
function get_ν_HH(E_bind::Real, mass::Real)::Real
    return √(2 * GRAIN_SURFACE_SITE_DENSITY * E_bind * K_B_SI / (mass * AMU_SI * π^2))
end

"""
    get_thermal_desorption_rate_constant(E_bind, ν_des, T_dust)

The thermal desorption rate constant,

``\\nu_{\\mathrm{des}}e^{-E_{\\mathrm{bind}}/T_{\\mathrm{dust}}}``
"""
function get_thermal_desorption_rate_constant(E_bind::Real, ν_des::Real, T_dust::Real)::Real
    return ν_des * get_arrhenius_probability(E_bind, T_dust)
end

"""
    get_thermal_diffusion_rate_constant(E_diff, ν_diff, T_dust)

The thermal diffusion rate constant,

``\\nu_{\\mathrm{diff}}e^{-E_{\\mathrm{diff}}/T_{\\mathrm{dust}}}``
"""
function get_thermal_diffusion_rate_constant(E_diff::Real, ν_diff::Real, T_dust::Real)::Real
    return ν_diff * get_arrhenius_probability(E_diff, T_dust)
end

function get_chemical_desorption_probability()
    error("unimplemented")
end

function get_enthalpy_of_reaction()
    error("unimplemented")
end

"""
    get_langmuir_hinshelwood_rate_constant(E_reac, μ_ij, T_dust, E_bind_i, E_bind_j, E_diff_i, E_diff_j, mass_i, mass_j, barrier_width=1.4)

Equation 21 from Vermarien et al, 2026.
"""
function get_langmuir_hinshelwood_rate_constant(
    E_reac::Real,
    T_dust::Real,
    E_bind_i::Real,
    E_bind_j::Real,
    E_diff_i::Real,
    E_diff_j::Real,
    mass_i::Real,
    mass_j::Real,
    μ_ij::Union{Real, Missing} = missing,
    barrier_width::Real = 1.4,
)::Real
    if ismissing(μ_ij)
        μ_ij = get_reduced_mass(mass_i, mass_j)
    end
    ν_HH_i = get_ν_HH(E_bind_i, mass_i)
    k_diff_i = get_thermal_diffusion_rate_constant(E_diff_i, ν_HH_i, T_dust)
    k_des_i = get_thermal_desorption_rate_constant(E_bind_i, ν_HH_i, T_dust)

    ν_HH_j = get_ν_HH(E_bind_j, mass_j)
    k_diff_j = get_thermal_diffusion_rate_constant(E_diff_j, ν_HH_j, T_dust)
    k_des_j = get_thermal_desorption_rate_constant(E_bind_j, ν_HH_j, T_dust)

    P = get_reaction_probability(E_reac, μ_ij, T_dust, barrier_width)
    κ_ij = max(ν_HH_i, ν_HH_j) * P

    f_comp = κ_ij / (κ_ij + k_diff_i + k_des_i + k_diff_j + k_des_j)

    # Should maybe be scaled some more?
    return f_comp * (k_diff_i + k_diff_j)
end

"""
    get_eley_rideal_rate_constant(freezeout_rate, γ, μ_ij, T_dust)

Equation 25 from Vermarien et al, 2026.
"""
function get_eley_rideal_rate_constant(
    freezeout_rate::Real,
    γ::Real,
    T_dust::Real,
    total_surface_abundance::Real,
    μ_ij::Real,
)::Real
    return freezeout_rate *
        get_reaction_probability(γ, μ_ij, T_dust) /
        total_surface_abundance
end

function get_crp_desorption_rate()
    error("unimplemented")
end
