@testset "Potential" begin
    # CarrierCapture.harmonic potential generation and finding the crossing
    Qmin = -25
    Qmax = 35
    NQ = 3000
    Q = range(Qmin, stop=Qmax, length=NQ)
    nev = 30
    hw = 0.02

    pot1 = potential(); pot1.name = "pot1"
    pot1.Q0 = 0; pot1.E0 = 1
    pot1.nev = nev
    pot1.func = x -> CarrierCapture.harmonic(x, hw; E₀ = pot1.E0, Q₀ = pot1.Q0)
    pot1.Q = Q
    pot1.E = pot1.func.(Q)

    pot2 = potential(); pot2.name = "pot2"
    pot2.Q0 = 2; pot2.E0 = 0
    pot2.nev = nev
    pot2.func = x -> CarrierCapture.harmonic(x, hw; E₀ = pot2.E0, Q₀ = pot2.Q0)
    pot2.Q = Q
    pot2.E = pot2.func.(Q)

    solve_pot!(pot1)
    solve_pot!(pot2)

    # eigenvalue test
    @test isapprox(pot1.ϵ, 1 .+ pot2.ϵ ; rtol=1E-5)

    Qx, Ex = find_crossing(pot1, pot2)
    @test Qx ≈ -4.22519859246400
    @test Ex ≈ 1.8541447195665985


    # Spline fit test
    nev = 60
    name = "D"

    Q1 = [30.66721918 29.860133 29.05306268 28.24612992 27.43911044 26.63197583 25.82507425 25.01797262 24.21096115 23.40387798 22.59690043 21.78991117 20.98281579 20.17570172 19.36884219 18.56169249 17.75463179 16.94772679 16.14061031 15.33347439 14.52663309 13.71945696 12.91240658 12.10544441 11.29841968 10.4913338 9.684370388 8.877289725 8.070184138]
    E1 = [8.0902 7.5970 7.0749 6.5242 5.9461 5.3451 4.7300 4.1147 3.5182 2.9637 2.4769 2.0819 1.7972 1.6315 1.5800 1.6237 1.7301 1.8586 1.9687 2.0283 2.0190 2.0673 1.9910 2.0528 1.9730 2.0857 2.4550 3.1653 4.3448]

    pot = potential(); pot.name = name
    pot.nev = nev
    pot.Q0 = Q1[findmin(E1)[2]]; pot.E0 = 1.69834
    pot.QE_data = DataFrame(Q = Q1[:], E = E1[:])
    pot.QE_data.E .+= - minimum(pot.QE_data.E) + pot.E0
    pot.Q = Q

    pot.func_type = "spline"
    ## spline fitting parameters
    params = Dict()
    params["order"] = 4
    params["smoothness"] = 0.001
    params["weight"] = [1 1 1 1 1 0.5 0.4 0.4 0.5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]

    pot.params = params

    fit_pot!(pot)
    @test pot.E == pot.func(Q)

    fit_pot!(pot)
    @test pot.E == pot.func(Q)

    pot.func_type = "bspline"
    fit_pot!(pot)
    @test pot.E == pot.func(Q)

    pot.func_type = "polyfunc"
    params["poly_order"] = 4
    params["p0"] = [1, 2, 3, 4, 6]
    fit_pot!(pot)
    @test pot.E ≈ pot.func(Q) atol = 1e-8

    pot.func_type = "harmonic_fittable"
    params["p0"] = [0.01]
    fit_pot!(pot)
    @test pot.E ≈ pot.func(Q) atol = 1e-8

    pot = potential_from_file("data/tio2_interpol.dat")
    @test pot.QE_data.E[2] ≈ -2336.08912022 atol = 1e-8
    # potential types
    # CarrierCapture.harmonic
    # double well
    # polyfunc
    # morse
    # get_spline
    # get_bspline

    #println(pot.T)

end
