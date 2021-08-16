import numpy as np
import shutil
import sys 
from params import *
from scipy import optimize
from scipy.fftpack import fft2,fftshift
from utils import *

def read_fields(include_relaxation):
    Emass_CBM = np.zeros((ngrid,ngrid))
    CBM = np.zeros((ngrid,ngrid))
    Emass_VBM = np.zeros((ngrid,ngrid))
    VBM = np.zeros((ngrid,ngrid))
    # reading the band edges from nscf.out
    for i in range(ngrid):
        for j in range(ngrid):
    
            directory = str(i) + '_' + str(j)
            #print(directory)
            PATH_Bilayer = os.path.join(ROOT_DIR, 'Bilayer', directory)
            os.chdir(PATH_Bilayer)
                    
            if os.path.exists(PATH_Bilayer + "/CRASH"):
               print('A calculation CRASHED in ' + PATH_Bilayer)         
    
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
                     
                        VBM[i,j] = 1000*energies[N_electrons - 1]
                        CBM[i,j] = 1000*energies[N_electrons]
                        
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
            
            klist = 2*np.asarray(klist) * np.pi * A_to_meters / alat
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
    
            Emass_VBM[i,j] = 2./(Emass_VB1 + Emass_VB2)
            Emass_CBM[i,j] = 2./(Emass_CB1 + Emass_CB2)
    
    VBM = np.array(symmetrize_c3(VBM))
    CBM = np.array(symmetrize_c3(CBM))
    Emass_VBM = np.array(symmetrize_c3(Emass_VBM))
    Emass_CBM = np.array(symmetrize_c3(Emass_CBM))
    
    plot_moire(VBM, os.path.join(ROOT_DIR, 'CM_Data', 'VBM_unrelax.pdf'),  'unit')
    plot_moire(CBM, os.path.join(ROOT_DIR, 'CM_Data', 'CBM_unrelax.pdf'),  'unit')
    plot_moire(Emass_VBM, os.path.join(ROOT_DIR, 'CM_Data', 'Emass_VBM_unrelax.pdf'),  'unit', 'emass')
    plot_moire(Emass_CBM, os.path.join(ROOT_DIR, 'CM_Data', 'Emass_CBM_unrelax.pdf'),  'unit', 'emass')
    
    VBMk = fft2(VBM)/ngrid**2
    CBMk = fft2(CBM)/ngrid**2
    Emass_VBMk = fft2(Emass_VBM)/ngrid**2
    Emass_CBMk = fft2(Emass_CBM)/ngrid**2
    
    plotterk(os.path.join(ROOT_DIR, 'CM_Data', 'VBM'), VBMk, 'GSFE', 'unrelax')
    plotterk(os.path.join(ROOT_DIR, 'CM_Data', 'CBM'), CBMk, 'GSFE', 'unrelax')
    plotterk(os.path.join(ROOT_DIR, 'CM_Data', 'Emass_inv_VBM'), Emass_VBMk, 'emass', 'unrelax')
    plotterk(os.path.join(ROOT_DIR, 'CM_Data', 'Emass_inv_CBM'), Emass_CBMk, 'emass', 'unrelax')
    
        
    if include_relaxation:
        ######## RELAXATION ###########
        print('Relaxation is included in the fields')
        u_rel = np.load(ROOT_DIR + '/Bilayer_relaxed.npy')
    
        VBM_relax, VBMk_relax = include_relax(VBMk, u_rel)
        CBM_relax, CBMk_relax = include_relax(CBMk, u_rel)
        Emass_VBM_relax, Emass_VBMk_relax = include_relax(Emass_VBMk, u_rel)
        Emass_CBM_relax, Emass_CBMk_relax = include_relax(Emass_CBMk, u_rel)
    
        moire_pot_VBM = np.max(VBM_relax)-np.min(VBM_relax)
        moire_pot_CBM = np.max(CBM_relax)-np.min(CBM_relax)        
        print('The moirè VBM potential is:', moire_pot_VBM)
        print('The moirè CBM potential is:', moire_pot_CBM)
                
        plot_moire(VBM_relax, os.path.join(ROOT_DIR, 'CM_Data', 'VBM_relax.pdf'),  'unit')
        plot_moire(CBM_relax, os.path.join(ROOT_DIR, 'CM_Data', 'CBM_relax.pdf'),  'unit')
        plot_moire(Emass_VBM_relax, os.path.join(ROOT_DIR, 'CM_Data', 'Emass_VBM_relax.pdf'),  'unit', 'emass')
        plot_moire(Emass_CBM_relax, os.path.join(ROOT_DIR, 'CM_Data', 'Emass_CBM_relax.pdf'),  'unit', 'emass')
    
        plotterk(os.path.join(ROOT_DIR, 'CM_Data', 'VBM'), VBMk_relax , 'GSFE', 'relax')
        plotterk(os.path.join(ROOT_DIR, 'CM_Data', 'CBM'), CBMk_relax, 'GSFE', 'relax')
        plotterk(os.path.join(ROOT_DIR, 'CM_Data', 'Emass_inv_VBM'), Emass_VBMk_relax, 'emass','relax')
        plotterk(os.path.join(ROOT_DIR, 'CM_Data', 'Emass_inv_CBM'), Emass_CBMk_relax, 'emass','relax')
    
    
