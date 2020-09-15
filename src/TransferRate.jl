include("Potential.jl")
export TransferCoord, get_coupling, get_reorg_energy, get_activation_energy

ħ = 6.582119514E-16 # epot⋅s
kB  = 8.6173303E-5 # epot⋅K⁻¹

"""
Stores three `Potential`s making up a transition between two diabatic states along an adiabatic potential energy surface.

## Fields

- `pot_d_i` and `pot_d_f` -- initial and final diabatic `Potential`s
- `pot_a` -- adiabatic `Potential`

## Constructor
    TransferCoord(pot_d_i::Potential, pot_d_i::Potential, pot_a::Potential)

"""
mutable struct TransferCoord
    # diabatic potentials
    pot_d_i::Potential; pot_d_f::Potential
    # adiabatic potential
    pot_a::Potential
end

"""
    get_coupling(tc::TransferCoord)

Calculate the electronic coupling between diabatic states

"""
function get_coupling(tc)
    crossing_Q, crossing_E = find_crossing(tc.pot_d_i, tc.pot_d_f)
    coupling = crossing_E - tc.pot_a.func(crossing_Q)

    return coupling
end

"""
    get_reorg_energy(tc::TransferCoord)

Calculate the reorganisation energy for the transfer

"""
function get_reorg_energy(tc::TransferCoord)
    min_i = tc.pot_d_i.E0
    min_f = tc.pot_d_f.E0
    # final diabatic potential evaulated at the initial minimum
    min_i_f = tc.pot_d_f.func(min_i)

    reorg_energy = min_i_f - min_f
    return reorg_energy
end

"""
    get_activation_energy(tc::TransferCoord)

Calculate the activation energy for the transfer

"""
function get_activation_energy(tc::TransferCoord)
    crossing_Q, crossing_E = find_crossing(tc.pot_d_i, tc.pot_d_f)
    activation_energy = crossing_E - tc.pot_d_i.E0

    return activation_energy
end

"""
    transfer_rate(tc::TransferCoord)

Calculate the transfer rate

"""
function transfer_rate(tc::TransferCoord)
    λ = get_reorg_energy(tc)
    Hab = get_coupling(tc)
    ΔG = get_activation_energy(tc)
    kBT = kB * tc.pot_a.T

    rate = (2*pi/ħ) * 1/sqrt(4*pi*λ*kBT) * Hab^2 * exp(-ΔG/kBT)

    return rate

end

