#!/usr/bin/env python3.7
# -*- coding=utf-8 -*-

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


def main(i_file, f_file):
    '''
    '''
    get_init_fin(i_file, f_file)



if __name__ == '__main__':
    '''
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i","--init",
                        help="initial input file (POSCAR format) ",default="./POSCAR_i")
    parser.add_argument("-f","--fin",
                        help="final input file (POSCAR format) ",default="./POSCAR_f")
    
    args = parser.parse_args()
    
    i_file = args.init
    f_file = args.fin

    main(i_file, f_file)

