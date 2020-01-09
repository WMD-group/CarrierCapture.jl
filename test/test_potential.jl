@testset "Potential" begin

    # Harmonic potential generation and finding the crossing 
    Qmin = -25
    Qmax = 35
    NQ = 3000
    Q = range(Qmin, stop=Qmax, length=NQ)
    nev = 30
    hw = 0.02

    pot1 = potential(); pot1.name = "pot1"
    pot1.Q0 = 0; pot1.E0 = 1
    pot1.nev = nev
    pot1.func = x -> harmonic(x, hw; E₀ = pot1.E0, Q₀ = pot1.Q0)
    pot1.Q = Q
    pot1.E = pot1.func.(Q)

    pot2 = potential(); pot2.name = "pot2"
    pot2.Q0 = 2; pot2.E0 = 0
    pot2.nev = nev
    pot2.func = x -> harmonic(x, hw; E₀ = pot2.E0, Q₀ = pot2.Q0)
    pot2.Q = Q
    pot2.E = pot2.func.(Q)
    
    solve_pot!(pot1)
    solve_pot!(pot2)

    # eigenvalue test
    @test isapprox(pot1.ϵ, 1 .+ pot2.ϵ ; rtol=1E-5)

    Qx, Ex = find_crossing(pot1, pot2)
    @test Qx ≈ -4.22519859246400
    @test Ex ≈ 1.8541447195665985

end
