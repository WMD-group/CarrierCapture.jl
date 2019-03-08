#!/usr/bin/env python3
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


def read_poscar(i_path, l_get_sorted_symbols=False):
    poscar = Poscar.from_file("{}".format(i_path))
    struct = poscar.structure
    if l_get_sorted_symbols:
        return struct, poscar.site_symbols
    else:
        return struct


def get_init_fin(i_file, f_file, disp_range=np.linspace(-1, 1, 11), output_dir='disp_dir'):
    '''
    '''
    # A. Alkauskas, Q. Yan, and C. G. Van de Walle, Physical Review B 90, 27 (2014)
    struct_i, sorted_symbols = read_poscar(i_file, True)
    struct_f, sorted_symbols = read_poscar(f_file, True)
    delta_R = struct_f.frac_coords - struct_i.frac_coords
    delta_R = (delta_R + 0.5) % 1 - 0.5
    
    lattice = struct_i.lattice.matrix #[None,:,:]
    delta_R = np.dot(delta_R, lattice)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Poscar(struct_i).write_file('disp_dir/POSCAR_i'.format(output_dir))
    # Poscar(struct_f).write_file('disp_dir/POSCAR_f'.format(output_dir))


    masses = np.array([spc.atomic_mass for spc in struct_i.species])
    delta_Q2 = masses[:,None] * delta_R ** 2

    print('Delta_Q^2', np.sqrt(delta_Q2.sum()))

    # print(delta_Q2.shape)
    # print(struct_i.lattice, struct_i.species, struct_i.cart_coords, )

    for frac in disp_range:
        disp = frac * delta_R
        struct = Structure(struct_i.lattice, struct_i.species, \
                           struct_i.cart_coords + disp, \
                           coords_are_cartesian=True)
        Poscar(struct).write_file('{0}/POSCAR_{1:03d}'.format(output_dir, int(np.rint(frac*10))))

def main(i_file, f_file, disp_range):
    '''
    '''
    get_init_fin(i_file, f_file, disp_range, output_dir='disp_dir_i')
    get_init_fin(f_file, i_file, disp_range, output_dir='disp_dir_f')



if __name__ == '__main__':
    '''
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i","--init",
                        help="initial input file (POSCAR format) ",default="./POSCAR_i")
    parser.add_argument("-f","--fin",
                        help="final input file (POSCAR format) ",default="./POSCAR_f")
    parser.add_argument("-d","--disp", nargs='+',
                        help="displacement range ",
                        default=[-1.0, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 0.6, 1.0])
    args = parser.parse_args()
    
    i_file = args.init
    f_file = args.fin
    disp_range = args.disp
 
    main(i_file, f_file, disp_range)

