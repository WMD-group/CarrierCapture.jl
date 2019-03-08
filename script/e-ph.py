import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Wavecar
np.set_printoptions(precision=4)


def read_wave_npy(wave_npy):
    return np.load(wave_npy)


def get_S_ifs(wf_path_i, wf_path_fs, verbose=False):
    """
    calculate overlaps
    real 
    Gamma point only but calculated by vasp_std
    we need a full list of C_nbk not reduced one.
    vasp_gam will return a half of plane waves and pymatgen will return Warning message.
    """
    S_if_list = []
    # initial a wave function
    wfc_i = read_wave_npy(wf_path_i) 
    wfc_i /= np.sqrt(np.vdot(wfc_i, (wfc_i)))

    for wf_path_f in wf_path_fs:
        # final wave functions(Q)
        wfc_f = read_wave_npy(wf_path_f)
        wfc_f /= np.sqrt(np.vdot(wfc_f, wfc_f))

        S_if = np.vdot(wfc_i, wfc_f)
        S_if_sq = (S_if*np.conjugate(S_if)).real
        S_if_list.append(np.sqrt(S_if_sq))

    print('S_if', S_if_list)

    return np.array(S_if_list)

def calc_W_if(ax, s_if, delta_eig, Qs, zero=-1):
    # calculate a dirivative of overlap betwenn initial end final wave functions
    # and plot and calculate e-ph coupling
    zero = zero if zero > -1 else int(len(f_files)/2)
    s_if = [s if i > zero else -s for i, s in enumerate(s_if) ]

    diff_overlap = np.gradient(s_if, Qs)

    ax.plot(Qs, s_if, '.-', label='S_if')
    ax.plot(Qs, diff_overlap, '.-', label='diff_overlap')
    print(delta_eig*diff_overlap[zero])

    ax.set_title('e capture')
    ax.legend()


def main(i_file, f_files, Q_max=7, delta_eig=1, zero=-1):
    delta_Q = Q_max/(len(f_files)) # assume uniform grid
    Qs = np.linspace(0, Q_max, len(f_files))
    
    zero = zero if zero > -1 else int(len(f_files)/2)
    
    s_if = get_S_ifs(i_file, f_files)
    s_if[zero] = 0
    fig, (ax1, ax2) = plt.subplots(1, 2)

    calc_W_if(ax1, s_if, delta_eig, Qs, zero)
    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--init",
                        help="initial wave function (npy format)")
    parser.add_argument("-f", "--fin", nargs='+',
                        help="final wave function s(npy format)")
    parser.add_argument("-d", "--deltaeig", type=float,
                        help="Delta eigenvalues of initial and final state")
    parser.add_argument("-z", "--zero", type=int,
                        help="index for Q_0")
    args = parser.parse_args()

    i_file = args.init
    f_files = args.fin
    zero = args.zero
    delta_eig = args.deltaeig
    
    main(i_file, f_files, delta_eig=delta_eig, zero=zero)

