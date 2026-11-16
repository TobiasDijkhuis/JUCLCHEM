using PhysicalConstants: AtomicMassConstants, PlankConstant, k_B, ħ

"""e^{-2a/ħ√{2μ_{ij}E_{reac}}}"""
function get_tunneling_probability(E_reac::Float64, barrier_width::Float64, μ_ij::Float64)::Float64
    return exp(-2 * barrier_width * 1.0e-10 / ħ * √(2 * μ_ij * AtomicMassConstants * E_reac * k_B))
end

"""e^{-E/T}"""
function get_arrhenius_probability(E_reac::Float64, temperature::Float64)::Float64
    return exp(-E_reac / temperature)
end
