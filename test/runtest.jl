# using NonradCapt    
include("../src/Potential.jl")

using Test

@testset "Vibrational Energy" begin
    etol = 1e-3
    NQ = 5000
    Qi = -5
    Qf = 5
    nev = 20
    dQ = (Qf-Qi)/NQ
    x = range(Qi, stop=Qf, length=NQ)
    ħω = 1 
    e, v = Potential.solve1D_ev_amu(x->Potential.harmonic(x, ħω), NQ=NQ , Qi=Qi, Qf=Qf, nev=nev)
    analytic_sol = UnitRange(0, nev) 
    @test isapprox(e, 0.5 .+  UnitRange(0, nev-1) ; rtol=etol)
end
