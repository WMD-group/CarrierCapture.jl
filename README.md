[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


CarrierCapture
==============

![](https://github.com/WMD-group/carriercapture/blob/master/schematics/Logo.png)

A set of codes to compute carrier capture and recombination rates in semiconducting compounds. 
This topic has a rich history starting from the work by [Huang and Rhys](http://rspa.royalsocietypublishing.org/content/204/1078/406.short) in 1950. 
Our implementation was inspired by the approach (and FORTRAN code) employed by [Alkauskas and coworkers](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.075202), but has been adapted
to also describe anharmonic potential energy surfaces. 

Installation
------------

The codes are written in [Julia](https://julialang.org), while the [Jupyter Notebooks](http://jupyter.org) also contain [Python](https://www.python.org), which are all assumed to be installed.

The [Brooglie](https://github.com/RedPointyJackson/Brooglie) package is used to solve the time-independent Schrödinger equation and can be installed by:

`Pkg.clone("https://github.com/RedPointyJackson/Brooglie")`  
`Pkg.add("Plots")`  
`Pkg.add("Arpack")`  
`Pkg.add("Polynomials")`  
`Pkg.add("CSV")`  

## Development

Development is in progress and hosted on [Github](https://github.com/WMD-group/carriercapture). 
Please use the [issue tracker](https://github.com/WMD-group/carriercapture/issues/) for feature requests, bug reports and more general questions. If you would like to contribute, please do so via a pull request.

Usage
-----

A typical usage will consist of three or four steps, implemented in a series of short programs:

1. Prepare a sequence of structures with displacements which interpolate between two defect states. Run single-point energy calculations on these structures, and extract the total energies.

2. Generate configuration coordinate diagrams with polynomial fits (`script here`). Solve these potential energy surfaces for the phonon wavefunctions for each defect and calculate the overlap between phonons in the donor and acceptor potentials.

3. Calculate the capture coefficient over a given temperature range

4. Calculate the lifetimes and rates for a given defect 

Code setup
----------
The capture.ipynb contains commands for each of the above steps, which call the Julia script `CaptureRate.jl`. 
    - general polynomial potential (calc_poly_wave_func)
    - calculating phonon overlap (calc_overlap)
    - calculating capture coefficient (calc_capt_coeff)

Input for calc_poly_wave_func: calc_poly_wave_func(potential_matrix_1, potential_matrix_2, poly_order, Qi=-10, Qf=10, NQ=100, nev=10, nev2=Nothing)
    - potential_matrix_1 & potential_matrix_2: potentials from files (named Potential1.txt and Potential2.txt).
    - Each potential is read from a file with Q and E in columns 1 and 2 respectively (energy in eV and Q in amu^(1/2) Angstrom). 
    - poly_order: order of the polynomial that you'd like to use to fit your potential data
    - Qi and Qf [amu^(1/2)*Å] define the domain where potential will be solved (Q ∈ [`Qi`,`Qf`]), discretised in `NQ` steps for each potential. 
    - `nev` the number of first energy levels of 1D potentials to solve. 
    - You can also include an additional argument nev2 = XX to solve for a different number of states in the second potential (by default nev2 = nev)

Input for calc_overlap: calc_overlap!(cc::CC; plt=Nothing, cut_off=0.25, σ=0.025, lplot=false)
    - cut_off: energetic difference criteria for overlap of phonons (Δϵ < cut_off)
    - σ: amount of smearing of delta functions for determining phonon overlap

Input for calc_capt_coeff: calc_capt_coeff(W, V, g, T_range, cc::CC)
    - V: volume of supercell [Å³]
    - g: configurational degeneracy 
    - W: electron-phonon coupling matrix element [ev/(amu^(1/2)*Å)]
    - T: temperature range for calculating capture coefficient (K)

Theory
------

The capture of electrons or holes by point defects in a crystalline materials requires the consideration of a number of factors, which include:

### Electronic matrix elements

> The electronic matrix element frequently causes feelings of discomfort (Stoneham, 1981)

Limitations
-----------

The expected use case for this code is to ... 

Examples
--------

The following examples are provided to illustrate some of the applications of these codes:

* GaAs defect (./example_GaAs) 


Citation
--------

If you use these codes in your own work, please cite...

Extended Reading List
---------------------

* [Heny and Lang, Nonradiative capture and recombination by multiphonon emission in GaAs and GaP (1977)](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.15.989) **Seminal contribution that introduces many concepts**

* [Huang, Adiabatic approximation theory and static coupling theory of nonradiative transition (1981)](http://engine.scichina.com/doi/10.1360/ya1981-24-1-27) **Useful context for the static approximation**

* [Stoneham, Non-radiative transitions in semiconductors (1981)](http://iopscience.iop.org/article/10.1088/0034-4885/44/12/001/meta) **Review on models and theory**

* [Markvart, Determination of potential surfaces from multiphonon transition rates (1981)](http://iopscience.iop.org/article/10.1088/0022-3719/14/15/002) **Treatment of anharmonicity**

* [Markvart, Semiclassical theory of non-radiative transitions (1981)](http://iopscience.iop.org/article/10.1088/0022-3719/14/29/006/meta) **Semiclassical treatment of matrix elements following Landau and Holstein**
