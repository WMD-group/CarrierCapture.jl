[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


# carriercapture

![](https://github.com/WMD-group/carriercapture/blob/master/schematics/Logo.png)

A set of codes to compute carrier capture and recombination rates in semiconducting compounds. 
This topic has a rich history starting from the work by [Huang and Rhys](http://rspa.royalsocietypublishing.org/content/204/1078/406.short) in 1950. 
Our implementation was inspired by the approach (and FORTRAN code) employed by [Alkauskas and coworkers](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.075202), but has been adapted
to also describe anharmonic potential energy surfaces. 

## Installation

The codes are written in [Julia](https://julialang.org), while the [Jupyter Notebooks](http://jupyter.org) also contain Python.
The [Brooglie](https://github.com/RedPointyJackson/Brooglie) package is used to solve the time-independent Schrödinger equation.

`Pkg.clone("https://github.com/RedPointyJackson/Brooglie")`

## Development

Development is in progress and hosted on [Github](https://github.com/WMD-group/carriercapture). 
Please use the [issue tracker](https://github.com/WMD-group/carriercapture/issues/) for feature requests, bug reports and more general questions. If you would like to contribute, please do so via a pull request.

## Theory

The capture of electrons or holes by point defects in a crystalline materials requires the consideration of a number of factors, which inlcude:

### Electronic matrix elements

> The electronic matrix element frequently causes feelings of discomfort (Stoneham, 1981)
