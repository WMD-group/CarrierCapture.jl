[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


# carriercapture

![](https://github.com/WMD-group/carriercapture/blob/master/schematics/Logo.png)

A set of codes to compute carrier capture and recombination rates in semiconducting compounds. 
This topic has a rich history starting from the work by [Huang and Rhys](http://rspa.royalsocietypublishing.org/content/204/1078/406.short) in 1950. 
Our implementation was inspired by the approach (and FORTRAN code) employed by [Alkauskas and coworkers](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.075202), but has been adapted
to also describe anharmonic potential energy surfaces. 

## Installation

The codes are written in [Julia](https://julialang.org), while the [Jupyter Notebooks](http://jupyter.org) also contain [Python](https://www.python.org), which are all assumed to be installed.

The [Brooglie](https://github.com/RedPointyJackson/Brooglie) package is used to solve the time-independent Schrödinger equation and can be installed by:

`Pkg.clone("https://github.com/RedPointyJackson/Brooglie")`
`Pkg.add("Plots")`
`Pkg.add("Arpack")`
`Pkg.add("Polynomials")`
`Pkg.add("DataFrames")`

## Development

Development is in progress and hosted on [Github](https://github.com/WMD-group/carriercapture). 
Please use the [issue tracker](https://github.com/WMD-group/carriercapture/issues/) for feature requests, bug reports and more general questions. If you would like to contribute, please do so via a pull request.

## Theory

The capture of electrons or holes by point defects in a crystalline materials requires the consideration of a number of factors, which inlcude:

### Electronic matrix elements

> The electronic matrix element frequently causes feelings of discomfort (Stoneham, 1981)

## Extended Reading List

* [Heny and Lang, Nonradiative capture and recombination by multiphonon emission in GaAs and GaP (1977)](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.15.989) **Seminal contribution that introduces many concepts**

* [Huang, Adiabatic approximation theory and static coupling theory of nonradiative transition (1981)](http://engine.scichina.com/doi/10.1360/ya1981-24-1-27) **Useful context for the static approximation**

* [Stoneham, Non-radiative transitions in semiconductors (1981)](http://iopscience.iop.org/article/10.1088/0034-4885/44/12/001/meta) **Review on models and theory**

* [Markvart, Determination of potential surfaces from multiphonon transition rates (1981)](http://iopscience.iop.org/article/10.1088/0022-3719/14/15/002) **Treatment of anharmonicity**

* [Markvart, Semiclassical theory of non-radiative transitions (1981)](http://iopscience.iop.org/article/10.1088/0022-3719/14/29/006/meta) **Semiclassical treatment of matrix elements following Landau and Holstein**

