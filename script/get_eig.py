#!/usr/bin/env python3
# -*- coding=utf-8 -*-

import sys
import os
import pickle

import argparse
import numpy as np

import xml.etree.cElementTree as ET
from monty.io import zopen, reverse_readfile

#from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin
from collections import defaultdict
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure

import matplotlib.pyplot as plt


E_MIN = -3
E_MAX = 3
dq = 0.01

# Copied from pymatgen
def _vasprun_float(f):
    """
    Large numbers are often represented as ********* in the vasprun.
    This function parses these values as np.nan
    """
    try:
        return float(f)
    except ValueError as e:
        f = f.strip()
        if f == '*' * len(f):
            warnings.warn('Float overflow (*******) encountered in vasprun')
            return np.nan
        raise e

# Copied from pymatgen
def _parse_varray(elem):
    if elem.get("type", None) == 'logical':
        m = [[True if i == 'T' else False for i in v.text.split()] for v in elem]
    else:
        m = [[_vasprun_float(i) for i in v.text.split()] for v in elem]
    return m

# Copied from pymatgen
def _parse_eigen(elem):
    eigenvalues = defaultdict(list)
    for s in elem.find("array").find("set").findall("set"):
        spin = Spin.up if s.attrib["comment"] == "spin 1" else Spin.down
        for ss in s.findall("set"):
            eigenvalues[spin].append(_parse_varray(ss))
    eigenvalues = {spin: np.array(v) for spin, v in eigenvalues.items()}
    elem.clear()
    return eigenvalues

def read_eigval(vasprun_path):
    filename='{}/vasprun.xml'.format(vasprun_path)
    #vr = Vasprun(filename)

    # to work with ISMEAR = -2
    for event, elem in ET.iterparse(zopen(filename, "rt")):
        if elem.tag == 'eigenvalues':
            eigval = _parse_eigen(elem)
            break
    return eigval

def plot_eigs(eigval, q):

    efermi = 0
    i_kpt = 0

    for spin, eigval_s in eigval.items():
        c = 'r' if float(spin) == 1 else 'b'
        for i_kpt in range(len(eigval_s)):
            eigval_sk = eigval[spin][i_kpt]
            occ = eigval_sk[:,1]
            y = eigval_sk[:,0]

            # occpied
            indx_occ = np.where(occ > 0.9)
            plt.plot([q+float(spin)*dq]*len(y[indx_occ]), y[indx_occ] - efermi, color=c, lw=0, marker='o')

            # unoccpied
            indx_unocc = np.where(occ < 0.3)
            plt.plot([q+float(spin)*dq]*len(y[indx_unocc]), y[indx_unocc] - efermi, color=c, lw=0, marker='o', fillstyle='none')

            # halfoccpied
            indx_hocc = (0.3 <= occ) & (occ <= 0.9)
            if len(indx_hocc) > 0:
                plt.plot([q+float(spin)*dq]*len(y[indx_hocc]), y[indx_hocc] - efermi, color=c, lw=0, marker='o', fillstyle='bottom')


def read_poscar(i_path, l_get_sorted_symbols=False):
    poscar = Poscar.from_file("{}".format(i_path))
    struct = poscar.structure
    if l_get_sorted_symbols:
        return struct, poscar.site_symbols
    else:
        return struct


def get_q(i_file, f_file, disp_range=None):
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
    return np.sqrt(delta_Q2.sum())

def plot(qs, eigvals, e_min, e_max):
    for q, eigval in zip(qs, eigvals):
        plot_eigs(eigval, q)
    plt.ylim((e_min, e_max))
    plt.show()

def read_data(paths, pivot_path):
    eigvals = []
    qs = []
    for path in paths:
        q = get_q('{}/POSCAR'.format(pivot_path), '{}/POSCAR'.format(path))
        eigval = read_eigval(path)
        eigvals.append(eigval)
        qs.append(q)
    return qs, eigvals

def save(qs, eigvals):
    with open('q_eig.pickle', 'wb') as handle:
        pickle.dump((qs, eigvals), handle, protocol=pickle.HIGHEST_PROTOCOL)

def main(paths, pivot_path, e_min, e_max):
    '''
    '''
    qs, eigvals = read_data(paths, pivot_path)
    save(qs, eigvals)
    plot(qs, eigvals, e_min, e_max)


if __name__ == '__main__':
    '''
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description="This script extracts and plots eigenvalues of a set of DFT (VASP) calculations")
    parser.add_argument("-p","--paths", nargs='+',
                        help="initial input file (POSCAR format) ",default=["DISP_-10", "DISP_-06", "DISP_-04", "DISP_-02", "DISP_-01", "DISP_000", "DISP_001", "DISP_002", "DISP_004", "DISP_006", "DISP_010"])
    parser.add_argument("-v","--pivot",
                        help="final input file (POSCAR format) ",default="DISP_-10")
    parser.add_argument("-e","--energy", nargs=2, type=float,
                        help="plot energy range ",default=[2, 5])

    args = parser.parse_args()

    paths = args.paths
    pivot_path = args.pivot
    e_min, e_max = args.energy

    main(paths, pivot_path, e_min, e_max)
