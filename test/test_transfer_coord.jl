@testset "Transfer rate" begin
    # adiabatic potential
    pot_adiabatic = pot_from_file("data/tio2_interpol_double.dat")
    pot_adiabatic.func_type = "spline"
    pot_adiabatic.Q0 = 0
    pot_adiabatic.E0 = 0

    fit_pot!(pot_adiabatic)

    # diabatic potentials
    pot_diabatic_i, pot_diabatic_f = cleave_pot(pot_adiabatic)
    pot_diabatic_i.func_type = "harmonic_fittable"
    pot_diabatic_f.func_type = "harmonic_fittable"
    pot_diabatic_i.E0 = pot_adiabatic.QE_data.E[1]
    pot_diabatic_i.Q0 = pot_adiabatic.QE_data.Q[1]
    pot_diabatic_f.E0 = pot_adiabatic.QE_data.E[end]
    pot_diabatic_f.Q0 = pot_adiabatic.QE_data.Q[end]

    filter_sample_points!(pot_diabatic_i)
    filter_sample_points!(pot_diabatic_f)

    fit_pot!(pot_diabatic_i)
    fit_pot!(pot_diabatic_f)

    tc = TransferCoord(pot_diabatic_i, pot_diabatic_f, pot_adiabatic)

    coupling = get_coupling(tc)
    @test coupling ≈ 0.050542 atol = 1e-6

    reorg_energy = get_reorg_energy(tc)
    @test reorg_energy ≈ 2.26961898898e6 atol = 1e-6

    activation_energy = get_activation_energy(tc)
    @test activation_energy ≈ 0.103882103 atol = 1e-6



end

