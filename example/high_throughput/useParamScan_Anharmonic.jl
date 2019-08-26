using Distributed
addprocs(length(Sys.cpu_info())-1) # use all cores
@everywhere using CarrierCapture
@everywhere using NPZ
@everywhere using SharedArrays

# params
@everywhere a_i = 11.7442
@everywhere a_f = 4.71261
@everywhere b_i = 0.0147757
@everywhere b_f = -0.0333741
@everywhere ΔQs = collect(range(0, stop=50, length=25)) # amu^0.5 Å
@everywhere ΔEs = collect(range(0, stop=2.5, length=10)) # eV. For solar materials, bandgap rarely exceeds 2eV

# initialise arrays
@everywhere ccArray = SharedArray{Float64}(length(ΔQs), length(ΔEs))

# Parallelise over the biggest loop
@sync @distributed for (qi,ΔQ) in collect(enumerate(ΔQs))
    for (ei,ΔE) in enumerate(ΔEs)
        poti, potf = fitMorseParams(a_i, a_f, b_i, b_f, ΔQ, ΔE)
        capt_coeff, barrier_height = getMorseCapture(poti, potf)

        for (i,cc) in enumerate(capt_coeff) # temperature dependence is stored in an array
            ccArray[qi, ei, i] = cc
        end

    end
end

npzwrite("ccArray.npz", ccArray)
