---
title: 'CarrierCapture.jl: Anharmonic Carrier Capture'
tags:
  - Julia
  - semiconductors
  - solar cells 
  - materials design
authors:
  - name: Sunghyun Kim
    orcid: 0000-0001-5072-6801
    affiliation: 1
  - name: Samantha N. Hood
    orcid: 0000-0001-9815-6289
    affiliation: 1  
  - name: Puck van Gerwen 
    affiliation: 1        
  - name: Lucy D. Whalley
    orcid: 0000-0002-2992-9871
    affiliation: 1      
  - name: Aron Walsh
    orcid: 0000-0001-5460-7033
    affiliation: "1, 2"
affiliations:
  - name: Department of Materials, Imperial College London, London, UK
    index: 1
  - name: Department of Materials Science and Engineering, Yonsei University, Seoul, Korea
    index: 2
date: 31 January 2020
bibliography: paper.bib
---
The modern theory of defects in semiconducting materials allows for the accurate prediction of equilibrium point defect concentrations and transitions levels within the band gap from first-principles quantum mechanical simulations [@park2018point]. The procedures for such calculations are now well established [@freysoldt2014first].

Beyond an equilibrium description, the operation and performance of optoelectronic devices, including solar cells and light-emitting diodes, relies on the interaction of points defects with non-equilibrium populations of electrons and holes. The non-radiative capture of charge carriers by lattice defects results in efficiency loss, which can range from supressing luminesence to irrereversible chemical degradation [@stoneham1981non].

`CarrierCapture.jl` is designed to calculate the rates of carrier capture by point defects from first-principles data. It builds on a large body of well established theory [@stoneham1981non], which was recently adapted to be compatible with quantities accessible from density functional theory (DFT) calculations [@alkauskas2014first]. In our implementation, we remove the harmonic approximation for the potential energy surfaces, which can be strongly asymmetric [@kim2019anham].

A standard workflow involves building potentials and computing carrier capture coefficients. `CarrierCapture.jl` provides handy tools: 

- Finding a best fit to the first-principles data.
- Solving one-dimensional Shrödinger equation for the potential energy surfaces.
- Computation of the overlap between nuclear wave functions.
- Computation of the capture coefficients as a function of temperature.  

We also provide auxiliary scripts to process the first-principles data.
To our best knowlege, `carrierCapture.jl` is the only open-source package that provide such functionalities.
Common computational workflows that reproduce published examples [@Kim2019kesterite; @kim2019anham] are available on Github pages. 
The API documentation including the guide to the intallation is also up-to-date on GitHub pages. 

**Interfacing to other codes:** A range of input paramaters are required to describe the bulk and defective properties of the material. Scripts are provided to extract these from DFT calculations using `VASP`, but these can easily be modified to read data from other packages. Caution should to be taken in checking covergence with respect to calculation settings, e.g. basis sets and k-point sampling, as small errors in relative energies can change the resulting carrier capture cross-sections by orders of magnitude. 


# Author contributions

[Sunghyun Kim](https://github.com/frssp) wrote the majority of the code base with contributions from [Samanth N. Hood](https://github.com/PaleBlueSam). [Lucy D. Whalley](https://github.com/lucydot) and [Puck van Gerwen](https://github.com/puckvg) performed detailed code testing and contributed to the example and test suite. All authors along with [Aron Walsh](https://github.com/aronwalsh) made decisions about code design and feature implementation. This manuscript was written with input from all co-authours.

# Acknowledgements

The development of this code has benefited through discussions with and contributions from many members of the Walsh research group including
Ji-Sang Park, Jarvist M. Frost, and Zijuan Xie, as well as John Buckeridge and Alexey A. Sokol from University College London. 

# References

