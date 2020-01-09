@testset "Vibrational Energy" begin
    NQ = 10000
    Qi = -4; Qf = 4
    nev = 10
    x = range(Qi, stop = Qf, length = NQ)
    ħω = 1 
    e, v = solve1D_ev_amu(x->harmonic(x, ħω,; E₀ = 0, Q₀ = 0), NQ = NQ , Qi = Qi, Qf = Qf, nev = nev)
    @test e ≈ 0.5 .+  UnitRange(0, nev-1) atol = 1e-3
end
