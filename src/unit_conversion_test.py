#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def harmonic(x, hbarw):
    # ev(amu)
    # ħ = 6.582119514E-16   # eV⋅s
    amu = 931.4940954E6   # eV / c^2
    hbarc = 0.19732697E-6 # eV m 
    a = amu / 2. * (hbarw/hbarc/1E10)**2
    return a*x*x

def hbarw2a(hbarw):
    amu = 931.4940954E6   # eV / c^2
    hbarc = 0.19732697E-6 # eV m 
    a = amu / 2. * (hbarw/hbarc/1E10)**2
    return a

def a2hbarw(a):
    amu = 931.4940954E6   # eV / c^2
    hbarc = 0.19732697E-6 # eV m 
    hbarw = (2 / amu * a)**0.5 *(hbarc*1E10) 
    return hbarw

def main():
    print(hbarw2a(0.04))
    print(a2hbarw(0.19138028580238858))


if __name__ == '__main__':
    main()