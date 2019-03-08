
# DX-
#  430        2.521943   1.000000
#  431        2.521943   1.000000
#  432        2.577513   1.000000
#  433        2.997799   1.000000
#  434        4.402615   0.000000
#  435        4.648891   0.000000


delta_eig=` bc <<< 4.402615-2.997799 `
python e-ph.py -i wf_CBM.npy -f wf_D_-02.npy wf_D_-01.npy wf_D_00{0..6}.npy -d ${delta_eig} -z 2 

delta_eig=` bc <<< -2.577513+2.997799 `
python e-ph.py -i wf_VBM.npy -f wf_D_-02.npy wf_D_-01.npy wf_D_00{0..6}.npy -d ${delta_eig} -z 2 
