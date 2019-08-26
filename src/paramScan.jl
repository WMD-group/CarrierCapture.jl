export getHarmonicQ_m, getMorseQ_m, fitHarmonicParams, fitMorseParams, getHarmonicCapture, getMorseCapture

function getHarmonicQ_m(ħω,E0)
    """
    This function gets the value of ΔQ corresponding to an activationless
    Marcus regime (Q_m) for two harmonic potential energy surfaces

    Parameters
    -----------
    ħω : energy of vibration of the initial state, in eV
    E0 : vertical shift between PES, in eV

    Returns
    -------
    Q_m : The value of ΔQ corresponding to an activationless Marcus regime, in amu^0.5 Å^-1
    """
    amu = 931.4940954E6   # eV / c^2
    ħc = 0.19732697E-6    # eV m

    a = amu/2 * (ħω/ħc/1E10)^2
    return (E0/a)^0.5
end

function getMorseQ_m(a, b, E0)
    """
    This function gets the value of Q corresponding to an activationless
    Marcus regime (Q_m) for two Morse potentials

    Parameters
    -----------
    a : first parameter of the Morse potential of the initial state
    b : second parameter of the Morse potential of the final state
    ΔE : vertical shift between PES, in eV

    Returns
    -------
    Q_m : The value of ΔQ corresponding to an activationless Marcus regime, in amu^0.5 Å^-1
    """

    Q_m = (1/b)*log(1-(E0/a)^0.5)
    return Q_m
end

function fitHarmonicParams(ħω_i, ħω_f, ΔQ, ΔE)
        # set nev
        nev_excited = 180
        nf_min = (1/(ħω_f))*((nev_excited+0.5)ħω_i+2ΔE)-0.5
        nev_relaxed = floor(nf_min+1)

        Q = range(-20-ΔQ, stop=20+ΔQ, length=5000) # ensure high enough to get good harmonic fit

        # potential 1
        potf = potential(); potf.name = "Relaxed state"
        potf.Q0 = ΔQ; potf.E0 = 0
        potf.nev = nev_relaxed # must overlap with initial excited state
        potf.func = x -> harmonic(x, ħω_f; E₀ = potf.E0, Q₀ = potf.Q0)
        potf.params["ħω"]= ħω_f
        potf.params["Q₀"]= ΔQ
	    potf.Q = Q
        potf.E = potf.func.(Q)
        solve_pot!(potf)

        # potential 2
        poti = potential(); poti.name = "Excited state"
        poti.Q0 = 0; poti.E0 = ΔE;
        poti.nev = nev_excited # must be converged for partition function
        poti.func = x -> harmonic(x, ħω_i; E₀ = poti.E0, Q₀ = poti.Q0)
        poti.Q = Q
        poti.E = poti.func.(Q)
        poti.params["ħω"]= ħω_i
        poti.params["ΔE"]= ΔE
	solve_pot!(poti)
        println("parameters fit to potentials")
        return poti, potf
end

function fitMorseParams(a_i, a_f, b_i, b_f, Q0, E0)
    """
    This function fits the parameters to the Morse potential energy surfaces and
    solves the potentials using a 1D S.E. solver

    Parameters
    ----------
    a_i: first parameter of the Morse potential of the initial state
    a_f: first parameter of the Morse potential of the final state
    b_i: second parameter of the Morse potential of the initial state
    b_f: seond parameter of the Morse potential of the final state
    Q0 : horizontal shift between the PES', in amu^0.5 Å^-1
    E0 : vertical shift between the PES', in eV

    Returns
    -------
    poti, potf: potential of initial and final states
    """
    # set number of eigenvalues (nev) to solve for in potential
    nev_excited = 180 # must be converged for partition function
    nev_relaxed = 800 # ensure max(ε_f) - max(ε_i) > E0

    Q = range(-20-Q0, stop=20+Q0, length=5000) # ensure high enough to get good harmonic fit

    # potential 1
    poti = potential(); poti.name = "Excited state"
    poti.Q0 = 0; poti.E0 = E0;
    poti.nev = nev_excited
    poti.func = x -> morse(x, [a_i, b_i]; E0 = poti.E0, Q0 = poti.Q0)
    poti.Q = Q
    poti.E = poti.func.(Q)
    poti.params["E0"]=E0
    poti.params["a"]= a_i
    poti.params["b"]= b_i
    solve_pot!(poti)

    # potential 2
    potf = potential(); potf.name = "Relaxed state"
    potf.Q0 = Q0; potf.E0 = 0
    potf.nev = nev_relaxed # must overlap with initial excited state
    potf.func = x -> morse(x, [a_f, b_f]; E0 = potf.E0, Q0 = potf.Q0)
    potf.params["a"]= a_f
    potf.params["b"]= b_f
    potf.params["Q0"]=Q0
    potf.Q = Q
    potf.E = potf.func.(Q)
    solve_pot!(potf)

    println("parameters fit to potentials")
    return poti, potf
end

function getHarmonicCapture(poti, potf)
    """
    Get the capture coefficient in a charge trapping process from an initial to
    final state described by harmonic potential energy surfaces.

    Parameters
    ----------
    poti: initial potential
    potf: final potential

    Returns
    -------
    capt_coeff: capture coefficient in units of cm^3 /s
    barrier_height: classical barrier to a capture process in eV
    """
    Volume = 1e-21 # cm^3
    σ = 0.8*potf.params["ħω"]
    g = 1 # degeneracy
    temperature = [300]

    ħω= potf.params["ħω"]
    E0 = poti.params["E0"]
    Q_m = getQ_m(ħω, E0)
    Q0 = potf.params["Q0"]
    W = 0.068/(Q0-Q_m) # e-ph coupling

    cc = conf_coord(poti, potf) # i, f
    println("conf coordinate fitted")
    cc.W = W
    cc.g = g
    calc_overlap!(cc; σ = σ) # assigns overlap matrix and δ matrix to cc
    println("overlap calculated")
    calc_capt_coeff!(cc, Volume, temperature)
    println("capture coefficient calculated")

    cross = try
                find_crossing(potf,poti)[2]
            catch
                50.0 # if it can't converge, set barrier very high
            end

    if cross==50.0
        barrier_height = 50.0
    else
        barrier_height = cross-poti.E0
    end
    capt_coeff = cc.capt_coeff[1]

    return capt_coeff, barrier_height
end

function getMorseCapture(poti, potf)
    """
    Get the capture coefficient in a charge trapping process from an initial to
    final state described by harmonic potential energy surfaces.

    Parameters
    ----------
    poti: initial potential
    potf: final potential

    Returns
    -------
    capt_coeff: capture coefficient in units of cm^3 /s
    barrier_height: classical barrier to a capture process in eV
    """
    Volume = 1e-21 # cm^3
    σ = 0.01
    g = 1 # degeneracy
    temperature = [300]

    Q0 = potf.params["Q0"]
    E0 = poti.params["E0"]
    a = potf.params["a"]
    b = potf.params["b"]
    Q_m = getMorseQ_m(a, b, E0)

    W = 0.068/(Q0-Q_m) # e-ph coupling

    cc = conf_coord(poti, potf) # i, f
    println("conf coordinate fitted")
    cc.W = W
    cc.g = g
    calc_overlap!(cc; σ = σ) # assigns overlap matrix and δ matrix to cc
    println("overlap calculated")
    calc_capt_coeff!(cc, Volume, temperature)
    println("capture coefficient calculated")

    cross = try
                find_crossing(potf,poti)[2]
            catch
                50.0 # if it can't converge, set barrier very high
            end

    if cross==50.0
        barrier_height = 50.0
    else
        barrier_height = cross-poti.E0
    end
    capt_coeff = cc.capt_coeff[1]

    return capt_coeff, barrier_height
end
