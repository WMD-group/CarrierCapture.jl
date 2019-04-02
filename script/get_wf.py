from pawpyseed.core.projector import * # also imports the wavefunction module
print('read wfs')
wf_02 = Wavefunction.from_directory('DISP_-02', setup_projectors=False)

wf_01 = Wavefunction.from_directory('DISP_-01', setup_projectors=False)
wf000 = Wavefunction.from_directory('DISP_000', setup_projectors=False)
wf001 = Wavefunction.from_directory('DISP_001', setup_projectors=False)
wf002 = Wavefunction.from_directory('DISP_002', setup_projectors=False)
wf003 = Wavefunction.from_directory('DISP_003', setup_projectors=False)
wf004 = Wavefunction.from_directory('DISP_004', setup_projectors=False)
print('realspace')

wf_CBM = wf000.write_state_realspace(433, 0, 0, dim=[336, 336, 336])
np.save('CBM.npy', wf_CBM)
print('CBM')
wf_VBM = wf000.write_state_realspace(431, 0, 0, dim=[336, 336, 336])
np.save('VBM.npy', wf_VBM)
print('VBM')

wf_defects = [ wf.write_state_realspace(432, 0, 0, dim=[336, 336, 336]) for wf in [wf_02, wf_01, wf000, wf001, wf002, wf003, wf004]]
print('DEFECTS')
np.save('D_-02.npy', wf_defects[0])
np.save('D_-01.npy', wf_defects[1])
np.save('D_000.npy', wf_defects[2])
np.save('D_001.npy', wf_defects[3])
np.save('D_002.npy', wf_defects[4])
np.save('D_003.npy', wf_defects[5])
np.save('D_004.npy', wf_defects[6])
