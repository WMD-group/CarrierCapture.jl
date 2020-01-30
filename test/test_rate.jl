@testset "Capture rate" begin

    Q = range(-20, stop=35, length=4250)

    # potential
    # pot1
    pot1 = potential(); pot1.name = "pot1"; pot1.nev = 30
    pot1.Q0 = 0; pot1.E0 = 1.5    
    pot1.func_type = "harmonic"
    params = Dict()
    params["hw"] = 0.03
    fit_pot!(pot1, Q; params = params)

    # pot2
    pot2 = potential(); pot2.name = "pot2"; pot2.nev = 50
    pot2.Q0 = 2; pot2.E0 = 0.6
    pot2.func_type = "harmonic"
    params = Dict()
    params["hw"] = 0.03
    fit_pot!(pot2, Q; params = params)

    solve_pot!(pot1)
    solve_pot!(pot2)

    # conf_coord
    Tmin = 10; Tmax = 800; NT = 100
    Volume = 1.28463E-21
    cut_off = 0.25; σ = 0.0075
    temperature = range(Tmin, stop=Tmax, length=NT)

    W = 0.205 # e-ph coupling
    g = 1 # degeneracy

    # build a configuration coordinate for the electron capture
    cc = conf_coord(pot1, pot2)
    cc.W = W; cc.g = g
    
    # calc_overlap!
    calc_overlap!(cc; cut_off = cut_off, σ = σ)

    using LinearAlgebra
    # zero-overlap between misaligned eigenvalues
    @test all(diag(cc.δ_matrix) .== 0) 

    # calc_capt_coeff!
    calc_capt_coeff!(cc, Volume, temperature)
    @test all(diag(cc.partial_capt_coeff[:, :, 1]) .== 0)
end
