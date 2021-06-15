import os
import numpy as np

ngrid = 3
n_atoms = 6
n_atoms_top = 3
alat = 3.152
geometry = 'hexagonal'
relax_input_file = 'qe_relax.in'
scf_input_file = 'qe_scf.in'
nscf_input_file = 'qe_nscf.in'
sbatch_file = 'qqe.job'
input_folder = 'input_files'
bandsx_file = 'qe_bandsx.in'

# effective mass related
nk_emass = 10
emass_frac = 0.05

K_edge = [0.0, 0.0, 0.0] #Location in K-space of the band edge in crystal coordinates 

def checkfolders():

    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

    if not os.path.exists(ROOT_DIR + '/Bilayer'):
           os.makedirs(ROOT_DIR + '/Bilayer')

    if not os.path.exists(ROOT_DIR + '/CM_Data'):
           os.makedirs(ROOT_DIR + '/CM_Data')
           
    return ROOT_DIR
    
    
def test_func(x, A0, A1):
    return A0 + A1 * (x**2)    
    
    
    
    
A_to_meters = 1e10
hbsquare = np.square(1.0545718*1e-34)
electron_mass = 9.11*1e-31
eV_to_J = 1.602177*1e-19    


