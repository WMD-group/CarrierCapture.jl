export TransferCoord, get_coupling, get_reorg_energy, get_activation_energy, get_transfer_rate, einstein_mobility

ħ = 6.582119514E-16 # eV⋅s
kB  = 8.6173303E-5 # eV⋅K⁻¹
ev_to_joule = 1.60218E-19
e_char = 1.60217662E-19 # C

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
    Q_min_i = tc.pot_d_i.Q0
    E_min_f = tc.pot_d_f.E0
    # final diabatic potential evaulated at the initial minimum
    E_min_i_f = tc.pot_d_f.func(Q_min_i)

    reorg_energy = E_min_i_f - E_min_f
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
function get_transfer_rate(tc::TransferCoord)
    λ = get_reorg_energy(tc)
    Hab = get_coupling(tc)
    ΔG = get_activation_energy(tc)
    kBT = kB * tc.pot_a.T

    rate = (2*pi/ħ) * 1/sqrt(4*pi*λ*kBT) * Hab^2 * exp(-ΔG/kBT)

    return rate
end

"""
    einstein_mobility(rate::Float64, n_neighbours::Int64, dist::Float64, temp::Float64)

## Fields

- `rate` -- Transfer rate
- `n_neighbours` -- Number of neighbours involved in the hopping
- `dist` -- Distance between sites in whichever unit
- `temp` -- Temperature in K

"""
function einstein_mobility(rate::Float64, n_neighbours::Int64, dist::Float64, temp::Float64)
    diffusion = dist^2 * n_neighbours * rate
    mob = e_char * diffusion / (kB * ev_to_joule * temp)
    return mob
end

