# Usage

A typical usage will consist of about three steps, implemented in a series of short programs which may be run from the command line. Input for the calculations is provided in `input.yaml`.

## 1. Preparation

Before `CarrierCapture`, you need to calculate potential energy surfaces of atomic vibrations (one-dimential Configuration Coordinate diagram; `1D-CC`) and _e-ph_ coupling matrix element (`W_if`). Prepare a sequence of structures with displacements which interpolate between two defect states. Run single-point energy calculations on these structures, and extract the total energies. Scripts for preprocessing can be found in `/script` which require a python library [`pymatgen`](http://pymatgen.org).

1. **Generate `1D-CC`**

   1. Calculate equilibirum geometries and total energies of defective supercells with charge states `q`(initial) and `q±1`(final) denoted `Conf.(q)` and `Conf.(q±1)`, respectively. 

   2. Generate interpolated and extrapolated structures between `Conf.(q)` (`POSCAR_i`) and `Conf.(q±1)` (`POSCAR_f`). You may use `gen_cc_struct.py`:

      ```bash
      $ gen_cc_struct.py -i POSCAR_i -f POSCAR_f -d -1 -0.6 -0.4 -0.2 -0.1 0 0.1 0.2 0.4 0.6 1.0 
      $ ls
      disp_dir_i
      disp_dir_f
      ```

   3. Run total-energy calculations for each structures. Example of the directory tree (`template` contains all input files for DFT calculations. Make sure DFT-program write wavefunctions (e.g. `LWAVE=.TRUE.` in `VASP`) for [the next stage `W_if`](#wif)):

      ```bash
      ├── 00_q2q±1
      │   ├── 00_q
      │   │   ├── DISP_000
      │   │   ├── DISP_001
      │   │   ├── ...
      │   │   ├── disp_dir -> ../../90_DISPLACEMENT/disp_dir_i
      │   │   └── template
      │   └── 01_q±1
      │       ├── DISP_000
      │       ├── DISP_001
      │       ├── ...
      │       ├── disp_dir -> ../../90_DISPLACEMENT/disp_dir_f
      │       └── template
      ├── 90_DISPLACEMENT
      │   ├── disp_dir_f
      │   │   ├── POSCAR_000
      │   │   ├── ...
      │   ├── disp_dir_i
      │   │   ├── POSCAR_000
      │   │   ├── ...
      ```

      You can submit jobs for all calculations using a following script in a high-performace computer with a batch system.

      ```bash
      #!/bin/bash -l
      
      for NUM in {-14,-13,-12,-11,-10,-09,-08,-07,-06,-05,-04,-03,-02,-01,000,001,002,003,004,005,006,007,008,009,010,011,012,013,014}
      do
         #if [ ! -d DISP_$NUM ]
         #then
           echo DISP_$NUM
           cp -r template DISP_$NUM
           cp disp_dir/POSCAR_$NUM DISP_$NUM/POSCAR
           cd  DISP_$NUM
           qsub run.pbs
           cd ../
         #fi
      done
      ```

   4. Calculate `𝛥Q` using `get_del_Q.py`. Generate `potential.csv` using following script.

      ```bash
      #!/bin/bash -l
      
      echo Q, E
      for NUM in {-14,-13,-12,-11,-10,-09,-08,-07,-06,-05,-04,-03,-02,-01,000,001,002,003,004,005,006,007,008,009,010,011,012,013,014}
      do
         en=`tail -n 1 DISP_$NUM/OSZICAR|awk '{print $3};'|sed  's/-./-0./' `
         del_Q=`get_del_Q.py -i $1 -f  DISP_$NUM/POSCAR`
         echo ${del_Q}, ${en}
      done
      ```

      Example `potential.csv` file:

      ```
      Q, E
      5.034147844442521, -0.28566871E+03
      4.027331082104984, -0.28581530E+03
      ...
      1.0068515567113572, -0.28581699E+03
      0.0, -0.28566641E+03
      ```

      

2. **Calcuate _e-ph_ coupling matrix element `W_if`** <a name="wif"></a>

   You already have eigenvalues, wave functions and configurations. Read [Work by Alkauskas and coworkers](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.075202) carefully.

   1. Find initial and final eigenvalues (`ϵ_i` and `ϵ_f`).

   2. Extract initial and final wavefunctions. For `VASP`, you can use `get_wf.py`:

      ```bash
      #!/bin/bash -l
      
      # Initial: VBM (hole capture)/CBM (electron capture)
      get_wf.py -f DISP_000/WAVECAR -k 0 -b <band_index_for_VBMorCBM> -s <spin_index(0or1)> -o wf
      
      mv wf_k0b<band_index>.npy wf_i.npy
      
      # Final: localized defect wavefucntion
      NBAND=<band_index_for_defect_wf>
      for NUM in {-10,-06,-04,-02,-01,000,001,002,004,006,010}
      do
          get_wf.py -f DISP_${NUM}/WAVECAR -k 0 -b ${NBAND} -s <spin_index(0or1)> -o wf_${NUM}
          mv wf_k0b<band_index>.npy wf_f_${NUM}.npy
      done
      ```

   3. Calculate the rate of change in overlap between initial and final wavefunctions `<ψ_i0|ψ_f(𝛥Q)>` as the geometry changes `𝛥Q`. The overlap can be calculated by using `get_overlap?.py`.

   4. ```W_if = (ϵ_f - ϵ_i) d<ψ_i0|ψ_f(𝛥Q)> / d𝛥Q```

