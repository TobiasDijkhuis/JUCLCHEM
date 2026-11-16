using PhysicalConstants.CODATA2022: AtomicMassConstant

const H2_FLUX_CT = 1.0e-10

"""
    get_thermal_velocity(mass, T)


``v_{\\mathrm{th}}=\\sqrt{\\frac{8k_\\mathrm{B}T}{πm}}``
"""
function get_thermal_velocity(mass::AbstractFloat, T::AbstractFloat)::AbstractFloat
    thermal_vel = √(8 * k_B * T / (π * AtomicMassConstant * mass))
    return thermal_vel
end

"""
    get_arrhenius_probability(E_reac, temperature)

``e^{-E/T}``"""
function get_arrhenius_probability(
    E_reac::AbstractFloat,
    temperature::AbstractFloat,
)::AbstractFloat
    return exp(-E_reac / temperature)
end

"""
    get_sticking_coefficient_Hollenbach_McKee(T_gas, T_dust)

Equation 3.7 from Hollenbach & McKee 1979.
"""
function get_sticking_coefficient_Hollenbach_McKee(
    T_gas::AbstractFloat,
    T_dust::AbstractFloat,
)::AbstractFloat
    S = 1 / (1 + 0.04 * √(T_gas + T_dust) + 0.2 * T_gas / 100 + 0.08 * (T_gas / 100)^2)
    return S
end

"""
    get_factor_1(μ, ν, E_H2, T_dust)

The first term in equation 15 from Cazaux & Tielens, 2002.
"""
function _get_factor_1(
    μ::AbstractFloat,
    ν::AbstractFloat,
    E_H2::AbstractFloat,
    T_dust::AbstractFloat,
)::AbstractFloat
    return μ * H2_FLUX_CT / (2 * ν * get_arrhenius_probability(E_H2, T_dust))
end

"""
    _get_factor_2(E_HC, E_S, E_HP, T_dust)

Equation 17 from Cazaux & Tielens, 2002."""
function _get_factor_2(
    E_HC::AbstractFloat,
    E_S::AbstractFloat,
    E_HP::AbstractFloat,
    T_dust::AbstractFloat,
)::AbstractFloat
    return (1 + √((E_HC - E_S) / (E_HP - E_S)))^2 / 4 *
           get_arrhenius_probability(E_S, T_dust)
end

"""
    _get_ϵ(\nu, E_HC, E_S, E_HP, T_dust

Equation 16 from Cazaux & Tielens, 2002.
"""
function _get_ϵ(
    ν::AbstractFloat,
    E_HC::AbstractFloat,
    E_S::AbstractFloat,
    E_HP::AbstractFloat,
    T_dust::AbstractFloat,
)::AbstractFloat
    return 1 / (
        1 +
        ν / (2 * H2_FLUX_CT) *
        get_arrhenius_probability(-1.5 * E_HC, T_dust) *
        (1 + √((E_HC - E_S) / (E_HP - E_S)))^2
    )
end

const GRAIN_DATA = ImmutableDict{
    "Si" => ImmutableDict{
        "μ" => 0.005,
        "E_S" => 110.0,
        "E_H2" => 320.0,
        "E_HP" => 450.0,
        "E_HC" => 3.0e4,
        "ν_H2" => 3.0e12,
        "ν_H" => 1.3e13,
        "CROSS_SECTION" => 8.437e22,
    },
    "Gr" => ImmutableDict{
        "μ" => 0.005,
        "E_S" => 260.0,
        "E_H2" => 520.0,
        "E_HP" => 800.0,
        "E_HC" => 3.0e4,
        "ν_H2" => 3.0e12,
        "ν_H" => 1.3e13,
        "CROSS_SECTION" => 0.5,
    },
}

function get_formation_efficiency(
    material::AbstractString,
    T_dust::AbstractFloat,
)::AbstractFloat
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

function get_CT_H2_form_efficiency(
    T_gas::AbstractFloat,
    T_dust::AbstractFloat,
)::AbstractFloat
    v_th = get_thermal_velocity(AtomicMassConstant, T_gas)
    sticking_coefficient = get_sticking_coefficient_Hollenbach_McKee(T_gas, T_dust)

    graphite_efficiency = get_formation_efficiency("Gr", T_dust)
    graphite_rate =
        v_th *
        sticking_coefficient *
        GRAIN_DATA["Gr"]["CROSS_SECTION"] *
        graphite_efficiency
    silicon_efficiency = get_formation_efficiency("Si", T_dust)
    silicon_rate =
        v_th * sticking_coefficient * GRAIN_DATA["Si"]["CROSS_SECTION"] * silicon_efficiency

    return 0.5 * (graphite_rate + silicon_rate)
end

"""
    get_tunneling_probability(E_reac, barrier_width, μ_ij)

``e^{-2a/\\hbar\\sqrt{2\\mu_{ij}E_{\\mathrm{reac}}}}``"""
function get_tunneling_probability(
    E_reac::AbstractFloat,
    barrier_width::AbstractFloat,
    μ_ij::AbstractFloat,
)::AbstractFloat
    return exp(
        -2 * barrier_width * 1.0e-10 / ħ * √(2 * μ_ij * AtomicMassConstant * E_reac * k_B),
    )
end

"""
    get_ν_HH(E_bind, mass)

``\\nu_{\\mathrm{HH}}=\\sqrt{2*SURFACE_SITE_DENSITY*E_{\\mathrm{bind}}/(π^2m)}``
"""
function get_ν_HH(E_bind::AbstractFloat, mass::AbstractFloat)::AbstractFloat
    return √(2 * SURFACE_SITE_DENSITY * E_bind * k_B / (mass * AtomicMassConstant * π^2))
end

"""
    get_thermal_desorption_rate_constant(E_bind, ν_des, T_dust)

The thermal desorption rate constant,

``\\nu_{\\mathrm{des}}e^{-E_{\\mathrm{bind}}/T_{\\mathrm{dust}}}``
"""
function get_thermal_desorption_rate_constant(
    E_bind::AbstractFloat,
    ν_des::AbstractFloat,
    T_dust::AbstractFloat,
)::AbstractFloat
    return ν_des * get_arrhenius_probability(E_bind, T_dust)
end

function get_chemical_desorption_probability()
    error("unimplemented")
end
