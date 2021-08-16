import os
import numpy as np
from utils import checkfolders

# Twist Angle
theta = 1.5
#Location in K-space of the band edge in crystal coordinates
K_edge = [0.0, 0.0, 0.0]  

ngrid = 9
n_atoms = 6
n_atoms_top = 3
alat = 3.192 # in A
C = 30.0 # in A, vacuum param
geometry = 'hexagonal'
relax_input_file = 'qe_relax.in'
relax_input_file_L1 = 'qe_relax_l1.in'
relax_input_file_L2 = 'qe_relax_l2.in'
scf_input_file = 'qe_scf.in'
nscf_input_file = 'qe_nscf.in'
sbatch_file = 'qqe.job'
input_folder = 'input_files'
bandsx_file = 'qe_bandsx.in'

# monolayer strain analysis
nstrain = 11
dstrain = 0.003

# effective mass related
nk_emass = 8
emass_frac = 0.05

#relaxation
kmax = 4.1 # expansion cutoff
q = 27 # upgrading real space fields

# CM related
Gcut = 3.1 # fields Fuorier components cutoff
eps = 1e-7
gdim = 4 # integer, CM energy cutoff

# bands
hspoints=[[1/3,1/3], [0.,0.], [0.5,0.], [1/3,-2/3]]
hsnames = ['K', '$\Gamma$', 'M', 'K']
Nk_points = 30
nbands = 12

#eigens/ldos
eigen_kpoints = {'G' : [0.,0.], 'K' : [1/3,1/3]}
nkmesh = 3


ROOT_DIR = checkfolders()

if geometry == 'hexagonal':
 
   R1 = alat * np.asarray([1., 0., 0.])
   R2 = alat * np.asarray([-1./2., np.sqrt(3.)/2., 0.])
   R3 = C * np.asarray([0.,0., 1.])
   
   lattice_matrix = np.vstack((R1,R2,R3))
   
   R1m = np.array([R1[1],-R1[0]])/np.radians(theta)
   R2m = np.array([R2[1],-R2[0]])/np.radians(theta)
   
   Scross = np.cross(R1[:2],R2[:2])
   G1 = 2*np.pi*np.array([R2[1],-R2[0]])/Scross 
   G2 = 2*np.pi*np.array([-R1[1],R1[0]])/Scross      
   
   G1m = G1 * np.radians(theta)
   G2m = G2 * np.radians(theta)

   uv1 = R1 / np.linalg.norm(R1[:2])
   uv2 = R2 / np.linalg.norm(R2[:2])
   
   Scross_norm = np.cross(uv1[:2],uv2[:2])   
   gv1 = 2*np.pi*np.array([uv2[1],-uv2[0]])/Scross_norm
   gv2 = 2*np.pi*np.array([-uv1[1],uv1[0]])/Scross_norm     

   uv1m = np.array([uv1[1],-uv1[0]])/np.radians(theta)
   uv2m = np.array([uv2[1],-uv2[0]])/np.radians(theta)
   gv1m = gv1* np.radians(theta)
   gv2m = gv2* np.radians(theta)
   MSL = np.matrix([uv1m,uv2m])

# various constants
A_to_bohr = 1.889725989
Ry_to_eV = 13.6056980659   
A_to_meters = 1e10
nm_to_meters = 1e9
hbsquare = np.square(1.0545718*1e-34)
electron_mass = 9.11*1e-31
eV_to_J = 1.602177*1e-19    

m0 =(0.5109989461*(10**9))/(8.987551787*(10**34))*1e-2 #in meV*s^2/A^2
hb2m0 = (6.582119514*(10**(-13)))**2/(2*m0) # !in meV


