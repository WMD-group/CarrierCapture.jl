@testset "Vibrational Energy" begin
    NQ = 10000; Qi = -4; Qf = 4; 
    Q = range(Qi, stop = Qf, length = NQ)

    ħω = 1
    func = x->CarrierCapture.harmonic(x, ħω,; E₀ = 0, Q₀ = 0)
    
    nev = 10
    e, v = CarrierCapture.solve1D_ev_amu(func, Q; nev = nev)
    
    @test e ≈ 0.5 .+  UnitRange(0, nev-1) atol = 1e-3
end
