using Distributed
addprocs(length(Sys.cpu_info())-1) # use all cores
@everywhere using CarrierCapture
@everywhere using NPZ
@everywhere using SharedArrays

# define the necessary parameters for the PES'
@everywhere ħω_i = 8e-3 # 8 meV
@everywhere ħω_f = 8e-3 # 8 meV
@everywhere ΔQs = collect(range(0, stop=25, length=25)) # amu^0.5 Å
@everywhere ΔEs = collect(range(0, stop=2.5, length=10)) # eV

# initialise array
@everywhere ccArray = SharedArray{Float64}(length(ΔQs), length(ΔEs))

# Parallelise over the biggest loop
@sync @distributed for (qi,ΔQ) in collect(enumerate(ΔQs))
    for (ei,ΔE) in enumerate(ΔEs)
        poti, potf = fitHarmonicParams(ħω_i, ħω_f, ΔQ, ΔE)
        capt_coeff, barrier_height = getHarmonicCapture(poti, potf)

        for (i,cc) in enumerate(capt_coeff) # temperature dependence is stored in an array
            ccArray[qi, ei, i] = cc
        end
    end
end

# store output in numpy array for analysis
npzwrite("ccArray.npz", ccArray)
