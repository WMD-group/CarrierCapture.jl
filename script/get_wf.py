import numpy as np
np.set_printoptions(4)

import argparse

def main(band_indice, defect_ind, dir_init, dir_final, kpoint_indice, spin_indice):
    from pawpyseed.core.projector import Wavefunction, Projector # also imports the wavefunction module
    wf_list = []
    overlap_list = []

    wf_basis = Wavefunction.from_directory(dir_init, setup_projectors=False)
    kpt_len = len(wf_basis.kpts)

    for dir_name in ["{}".format(_) for _ in dir_final]:
        wf_trap = Wavefunction.from_directory(dir_name, setup_projectors=False)
        pr = Projector(wf_trap, wf_basis)
        result = pr.single_band_projection(band_num=defect_ind)

        result = result.reshape((int(len(result)/kpt_len), kpt_len))
        result = result[:, int(kpoint_indice)]
        for b in band_indice:
            band_in=str(b)
        band_i = int(band_in)*int(spin_indice+1)
        overlap_list.append(result[band_i])

    print("\n\n\n")
    print("=================================")
    print("----------- Overlaps ------------")
    print("=================================")
    print(np.abs(np.array(overlap_list)))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Pawpyseed wrapper: calculating overlap between wave functions \
        e.g. python get_wf.py   -d 288 -b 291  -D DISP_000 -i DISP_-02 DISP_-01 DISP_001 DISP_002 -k 0 -s 1')

    parser.add_argument('-D','--dir_init', help='directory for the initial state', type=str, required=True, default="DISP_000")
    parser.add_argument('-i','--input', help='defect for deformed structures', type=str, nargs='+', required=True)

    parser.add_argument('-b','--band', help='band index for a initial state (start at 0)', nargs='+', type=int, required=True)
    parser.add_argument('-d','--defect', help='band index for a final state (start at 0)', type=int, required=True)
    parser.add_argument('-k','--kpoint',help='kpoint index (start at 0)', type=int, required=True, default=0)
    parser.add_argument('-s','--spin',help='spin-unpolarised/polarised (0/1) calculation', type=int, required=True, default=0)
    args = parser.parse_args()


    band_indice = args.band
    defect_ind = args.defect
    dir_init  = args.dir_init
    dir_final = args.input
    kpoint_indice = args.kpoint
    spin_indice = args.spin

    main(band_indice, defect_ind, dir_init, dir_final, kpoint_indice, spin_indice)
