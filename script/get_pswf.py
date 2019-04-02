#!/usr/bin/env python3
# -*- coding=utf-8 -*-

import sys
import argparse
import numpy as np
from pymatgen.io.vasp.outputs import Wavecar

def read_WAVECAR(file='WAVECAR', ks=0, bs=1, s=0, output='wf'):
    # The list of coefficients for each k-point and band for reconstructing the wavefunction. 
    # For non-spin-polarized, the first index corresponds to the kpoint and the second corresponds to the band 
    # (e.g. self.coeffs[kp][b] corresponds to k-point kp and band b). 
    # For spin-polarized calculations, the first index is for the spin.
    wavecar = Wavecar(file)
    for k in ks:
        for b in bs:
            #wf = wavecar.fft_mesh(kpoint=k, band=b-1)
            coeffs = np.array(wavecar.coeffs)
            if coeffs.ndim == 3:
                wf = coeffs[s][k][b]
            else:
                wf = coeffs[k][b]
            print('wf.shape', wf.shape)
            print('wavecar.kpoints', wavecar.kpoints)
            np.save('{}_k{}b{}.npy'.format(output, k, b), wf)


def main(file, s, k, b, output):
    read_WAVECAR(file, k, b, s, output)


if __name__ == '__main__':
    '''
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file",
                        help="VASP WAVECAR ", default="./WAVECAR")
    parser.add_argument("-o","--output", 
                    help="save ", default="wf")
    parser.add_argument("-b","--band", type=int, nargs='+',
                        help="band index starting from 1", default=[1])
    parser.add_argument("-k","--kpoint", type=int, nargs='+',
                        help="kpoint index starting from 0", default=[0])
    parser.add_argument("-s","--spin", default=0, type=int,
                        help="spin index 1 or 0",)

    args = parser.parse_args()
    file = args.file
    output = args.output
    band = args.band
    kpoint = args.kpoint
    spin = args.spin
    main(file, spin, kpoint, band, output)

