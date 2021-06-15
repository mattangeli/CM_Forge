import numpy as np
import shutil
import sys 
from params import *
from scipy import optimize
import matplotlib.pyplot as plt

ROOT_DIR = checkfolders()      

moire_potential_VB = []
moire_potential_CB = []
emass_VB_list = []
emass_CB_list = []

# reading the band edges from nscf.out
for i in range(ngrid):
    for j in range(ngrid):

        directory = str(i) + str(j)
        PATH_Bilayer = os.path.join(ROOT_DIR, 'Bilayer', directory)
        os.chdir(PATH_Bilayer)

     #   if os.path.exists(PATH_Bilayer + "/CRASH"):
     #      sys.exit('A calculation CRASHED in ' + PATH_Bilayer)         

# reading the band edge
        energies = []
        with open('nscf_edge.out', 'r') as infile:
             lines_edge = infile.readlines()
             
             for l,line in enumerate(lines_edge):

                 strings = line.split()
                 
                 if 'Cholesky' in strings:
                    sys.exit('Problem computing Cholesky decomposition in ' + PATH_Bilayer)   
              
                 if 'Kohn-Sham' in strings and 'states=' in strings:
                    N_states = int(strings[4])
              
                 if 'number' in strings and 'electrons' in strings:
                    N_electrons = int(float(strings[4]))  

                 if 'End' in strings and 'band' in strings and 'structure' in strings:    
                    if N_states % 8 == 0:
                       nrows = int(N_states / 8)
                    else:
                       nrows = int(N_states / 8) + 1
                       
                    for k in range(4, 4 + nrows):
                        estr = lines_edge[l+k].split()  
                        for e, eigen in enumerate(estr):
                            energies.append(float(eigen))                             

                    #sanity check
                    if len(energies) != N_states:
                       print('ERROR: wrong number of Kohn-Sham states')           
                 
                    moire_potential_VB.append(energies[N_electrons - 1])
                    moire_potential_CB.append(energies[N_electrons])
                   
                    break
        
        VB1 = []
        CB1 = []
        klist = []
        with open('emass_band1.gnu', 'r') as infile:
             lines_emass = infile.readlines()
             
             VB_idx = (nk_emass + 2) * (N_electrons - 1) + 1
             CB_idx = (nk_emass + 2) * (N_electrons) + 1

             for l in range(VB_idx, VB_idx + nk_emass):
                 strings = lines_emass[l].split()
                 
                 VB1.append(float(strings[1]))    
                 klist.append(float(strings[0]))
                    
             for l in range(CB_idx, CB_idx + nk_emass):
                 strings = lines_emass[l].split()
                 
                 CB1.append(float(strings[1]))            
        
        klist = 2*np.asarray(klist) *np.pi * A_to_meters / alat
        VB1 = np.asarray(VB1)*eV_to_J
        CB1 = np.asarray(CB1)*eV_to_J        
        
        params, params_covariance = optimize.curve_fit(test_func, klist, VB1, p0=[0., -1e-38])  
        Emass_VB1 = (hbsquare/(2*params[1] * electron_mass))

        params, params_covariance = optimize.curve_fit(test_func, klist, CB1, p0=[0., 1e-38])  
        Emass_CB1 = (hbsquare/(2*params[1] * electron_mass))
        
        VB2 = []
        CB2 = []
        klist = []
        with open('emass_band2.gnu', 'r') as infile:
             lines_emass = infile.readlines()
             
             VB_idx = (nk_emass + 2) * (N_electrons - 1) + 1
             CB_idx = (nk_emass + 2) * (N_electrons) + 1

             for l in range(VB_idx, VB_idx + nk_emass):
                 strings = lines_emass[l].split()
                 
                 VB2.append(float(strings[1]))    
                 klist.append(float(strings[0]))
                    
             for l in range(CB_idx, CB_idx + nk_emass):
                 strings = lines_emass[l].split()
                 
                 CB2.append(float(strings[1]))    
        
        os.chdir(ROOT_DIR)
        
        klist = 2*np.asarray(klist) * np.pi * A_to_meters / alat
        VB2 = np.asarray(VB2) * eV_to_J
        CB2 = np.asarray(CB2) * eV_to_J        
        
        params, params_covariance = optimize.curve_fit(test_func, klist, VB1, p0=[0., -1e-38])  
        Emass_VB2 = (hbsquare/(2*params[1] * electron_mass))

        params, params_covariance = optimize.curve_fit(test_func, klist, CB1, p0=[0., 1e-38])  
        Emass_CB2 = (hbsquare/(2*params[1] * electron_mass))        

        emass_VB_list.append((Emass_VB1 + Emass_VB2)/2.)
        emass_CB_list.append((Emass_CB1 + Emass_CB2)/2.)

GSFE = np.vstack((moire_potential_VB, moire_potential_CB)).T
Emass = np.vstack((emass_VB_list, emass_CB_list)).T
PATH_GSFE = os.path.join(ROOT_DIR, 'CM_Data', 'GSFE_unrelax.dat')
PATH_Emass = os.path.join(ROOT_DIR, 'CM_Data', 'Emass_unrelax.dat')
np.savetxt(PATH_GSFE, GSFE)
np.savetxt(PATH_Emass, Emass)



