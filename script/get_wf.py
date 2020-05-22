import numpy as np
np.set_printoptions(4)

import argparse

def main(band_indice, defect_ind, dir_init, dir_final):
    from pawpyseed.core.projector import Wavefunction, Projector # also imports the wavefunction module
    wf_list = []
    overlap_list = []

    wf_basis = Wavefunction.from_directory(dir_init, setup_projectors=False)
    kpt_len = len(wf_basis.kpts)

    for dir_name in ["{}".format(_) for _ in dir_final]:
        wf_trap = Wavefunction.from_directory(dir_name, setup_projectors=False)
        pr = Projector(wf_trap, wf_basis)
        result = pr.single_band_projection(band_num=defect_ind)

        # result = result.reshape((kpt_len, int(len(result)/kpt_len)))
        result = result.reshape((int(len(result)/kpt_len), kpt_len))
        result = result[:, 0]
        overlap_list.append(result[band_indice])

    print("\n\n\n")
    print("=================================")
    print("----------- Overlaps ------------")
    print("=================================")
    print(np.abs(np.array(overlap_list)))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Pawpyseed wrapper: calculating overlap between wave functions \
        e.g. python get_wf.py   -d 288 -b 291  -D DISP_000 -i DISP_-02 DISP_-01 DISP_001 DISP_002 ')

    parser.add_argument('-D','--dir_init', help='directory for the initial state', type=str, required=True, default="DISP_000")
    parser.add_argument('-i','--input', help='defect for deformed structures', type=str, nargs='+', required=True)

    parser.add_argument('-b','--band', help='band index for a initial state', nargs='+', type=int, required=True)
    parser.add_argument('-d','--defect', help='band index for a final state', type=int, required=True)
    args = parser.parse_args()


    band_indice = args.band
    defect_ind = args.defect
    dir_init  = args.dir_init
    dir_final = args.input

    main(band_indice, defect_ind, dir_init, dir_final)
