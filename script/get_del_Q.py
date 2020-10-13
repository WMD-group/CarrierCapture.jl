#!/usr/bin/env python3.7
# -*- coding=utf-8 -*-
"""
Get the mass weighted coordinate difference between two geometries.

The distance can be in a straight line, or projected on to a vector of initial and final geometries.

"""

import sys
import os

import argparse
import numpy as np

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Poscar, Kpoints
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure

import matplotlib.pyplot as plt

def sort_structure(structure, specie_order):
    ''' species sort
    '''
    d = dict([(ele, i)for i, ele in enumerate(specie_order)])

    return structure.get_sorted_structure(key=lambda x:d[x.species_string])


def read_poscar(i_path, l_get_sorted_symbols=False):
    poscar = Poscar.from_file("{}".format(i_path))
    struct = poscar.structure
    if l_get_sorted_symbols:
        return struct, poscar.site_symbols
    else:
        return struct


def get_init_fin(i_file, f_file, disp_range=None):
    '''
    '''
    # A. Alkauskas, Q. Yan, and C. G. Van de Walle, Physical Review B 90, 27 (2014)


    struct_i, sorted_symbols = read_poscar(i_file, True)
    struct_f, sorted_symbols = read_poscar(f_file, True)
    delta_R = struct_f.frac_coords - struct_i.frac_coords
    delta_R = (delta_R + 0.5) % 1 - 0.5
    
    lattice = struct_i.lattice.matrix
    delta_R = np.dot(delta_R, lattice)

    # delta_R = struct_i.cart_coords - struct_f.cart_coords

    masses = np.array([spc.atomic_mass for spc in struct_i.species])
    delta_Q2 = masses[:,None] * delta_R ** 2
    print(np.sqrt(delta_Q2.sum()))
    # delta_Q = np.sqrt(masses[:,None])* delta_R 
    # print(delta_Q.sum())


def main(args):
    struct_i = read_poscar(args.init)
    struct_f, sorted_symbols = read_poscar(args.fin, True)
    delta_R = struct_f.frac_coords - struct_i.frac_coords
    delta_R = (delta_R + 0.5) % 1 - 0.5
    lattice = struct_i.lattice.matrix
    delta_R = np.dot(delta_R, lattice)

    masses = np.array([spc.atomic_mass for spc in struct_i.species])

    # allow for no mass weighting
    if args.no_weight:
        delta_Q2 = delta_R ** 2
    else:
        delta_Q2 = masses[:,None] * delta_R ** 2
    delta_Q = np.sqrt(delta_Q2.sum())

    print("Delta Q:", delta_Q)

    if args.med != "unassigned":
        struct_m = read_poscar(args.med)
        delta_M = struct_m.frac_coords - struct_i.frac_coords
        delta_M = (delta_M + 0.5) % 1 - 0.5
        delta_M = np.dot(delta_M, lattice)
        # project the midpoint on the delta_R vector
        # einsum for row-wise dot product
        delta_M_proj = np.einsum('ij,ij->i',delta_R, delta_M)/np.linalg.norm(delta_R, axis=1)
        if args.no_weight:
            delta_M_proj_Q2 = delta_M_proj ** 2
        else:
            #delta_M_proj_Q2 = masses[:,None] * delta_M_proj ** 2
            delta_M_proj_Q2 = masses * delta_M_proj ** 2
        delta_M_proj_Q = np.sqrt(delta_M_proj_Q2.sum())

        print("Projected delta Q:", delta_M_proj_Q)
        print("Fractional displacement:", delta_M_proj_Q/delta_Q)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description="This script calculates the atomic mass weighted distance between two structures.")
    parser.add_argument("-i","--init",
                        help="initial input file (POSCAR format) ",default="./POSCAR_i")
    parser.add_argument("-f","--fin",
                        help="final input file (POSCAR format) ",default="./POSCAR_f")
    parser.add_argument("-m","--med",nargs="?",default='unassigned',
            help="optional medium geometry (POSCAR format) to get its the fractional displacement between the initial and final geometries")
    parser.add_argument("-nw","--no_weight",action="store_true",help="Turn off mass weighting")

    args = parser.parse_args()

    main(args)

