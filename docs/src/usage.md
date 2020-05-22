# Usage

A typical usage will consist of about three steps; 1. preparation, 2. building `potential`, and 3. computing capture coefficient. Find more detail for step 2 and step 3 in [example notebooks](https://github.com/WMD-group/CarrierCapture.jl/tree/master/example/notebook). The command line interface is depreciated and not recommended.

## 1. Preparation

Before `CarrierCapture`, you need to calculate potential energy surfaces of atomic vibrations (one-dimensional Configuration Coordinate diagram; `1D-CC`) and _e-ph_ coupling matrix element (`W_if`). Prepare a sequence of structures with displacements which interpolate between two defect states. Run single-point energy calculations on these structures, and extract the total energies. Scripts for preprocessing can be found in `/script` which require the [`pymatgen`](http://pymatgen.org) python library.

1. **Generate `1D-CC`**

   1. Calculate equilibrium geometries and total energies of defective supercells with charge states `q`(initial) and `q±1`(final) denoted `Conf.(q)` and `Conf.(q±1)`, respectively.

   2. Generate interpolated and extrapolated structures between `Conf.(q)` (`POSCAR_i`) and `Conf.(q±1)` (`POSCAR_f`). You may use `gen_cc_struct.py`:

      ```bash
      $ gen_cc_struct.py -i POSCAR_i -f POSCAR_f -d -1 -0.6 -0.4 -0.2 -0.1 0 0.1 0.2 0.4 0.6 1.0
      $ ls
      disp_dir_i
      disp_dir_f
      ```

   3. Run total-energy calculations for each structure. Example of the directory tree (`template`) contains all input files for DFT calculations. Make sure DFT-program writes wavefunctions (e.g. `LWAVE=.TRUE.` in `VASP`) for [the next stage `W_if`](#wif)):

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

      You can submit jobs for all calculations using a following script in a high-performance computer with a batch system.

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

   4. Calculate `𝛥Q` using `get_del_Q.py`. Generate `potential.csv` using following script.<a name="qe_data"></a>

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

   You already have eigenvalues, wavefunctions and configurations. Read [Work by Alkauskas and coworkers](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.075202) carefully.

   1. Find initial and final eigenvalues (`ϵ_i` and `ϵ_f`).

   2. Calculate overlap initial and final wavefunctions. For `VASP`, you can use `get_wf.py`:
   
      ```bash
      $ python get_wf.py   -d 288 -b 291  -D DISP_000 -i DISP_-02 DISP_-01 DISP_001 DISP_002
      
      GRID ENCUT 918.2911873931497
      finished making projector list
      --------------
      ran get_projector_list in 0.033210 seconds
      ---------------
      STARTING PROJSETUP
      started setup_proj
      calculating projector_values
      onto_projector calcs
      Done
      --------------
      
      ...
      
      =================================
      ----------- Overlaps ------------
      =================================
      [[0.0574]
       [0.0267]
       [0.0345]
       [0.0643]]
      ```
   
   3. Calculate the rate of change in overlap between initial and final wavefunctions `<ψ_i0|ψ_f(𝛥Q)>` as the geometry changes `𝛥Q`:
      `W_if = (ϵ_f - ϵ_i) d<ψ_i0|ψ_f(𝛥Q)> / d𝛥Q`. See more detail in [this example](https://github.com/WMD-group/CarrierCapture.jl/blob/master/example/notebook/e-ph.ipynb).
   
## 2. Building `potential`
See [Example](https://github.com/WMD-group/CarrierCapture.jl/blob/master/example/notebook/Anharmonic%20(DX%20center).ipynb).

  1. Use `fit_pot!` to find a best fit to the [data of `Q` and `E`](#qe_data). 
  2. Use `solve_pot!` to solve 1D Shrödinger equation for the potential energy surface (PES).

## 3. Computing capture rates
See [Example](https://github.com/WMD-group/CarrierCapture.jl/blob/master/example/notebook/Anharmonic%20(DX%20center).ipynb).

  1. Use `calc_overlap!` to calculate the overlap between phonon wave functions.
  2. Use `calc_capt_coeff!` to calculate the capture coefficient as a function of temperature.



# High-Throughput Usage
High-throughput usage is possible by preparing files in a similar method to the examples, `useParamScan_Harmonic.jl` and `useParamScan_Anharmonic.jl`. It is recommended that a high-performance computer rather than a personal machine is used, depending on how many calculations are performed. The code can then be run remotely using `nohup julia useParamScan_Harmonic &`.

The steps are as follows,
1. Define the parameters for which capture coefficient C will be calculated.

2. Calculate C for these parameters, in parallel over the largest parameter range (usually ΔQ)

3. Find `ccArray.npz` to analyse the capture coefficient as a function of its parameters.

