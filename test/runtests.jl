  
using CarrierCapture
using Test

@testset "Vibrational Energy" begin
    etol = 1e-3
    NQ = 5000
    Qi = -5
    Qf = 5
    nev = 20
    dQ = (Qf-Qi)/NQ
    x = range(Qi, stop = Qf, length = NQ)
    ħω = 1 
    e, v = solve1D_ev_amu(x->harmonic(x, ħω,; E₀ = 0, Q₀ = 0), NQ = NQ , Qi = Qi, Qf = Qf, nev = nev)
    @test isapprox(e, 0.5 .+  UnitRange(0, nev-1) ; rtol=etol)
end
