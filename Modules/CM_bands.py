import numpy as np
import shutil
import sys
import os 
from params import *
from utils import *
from itertools import product

def continuum_model_bands(include_emass_correction, plot_eigen, plot_ldos):   
   
   print('computing the bandstructure for the twist angle: ', theta)
   ROOT_DIR = checkfolders()      
   PATH_VBM = os.path.join(ROOT_DIR, 'CM_Data', 'VBM_FT.dat')
   PATH_CBM = os.path.join(ROOT_DIR, 'CM_Data', 'CBM_FT.dat')
   PATH_EMASS_VBM = os.path.join(ROOT_DIR, 'CM_Data', 'Emass_inv_VBM_FT.dat')
   PATH_EMASS_CBM = os.path.join(ROOT_DIR, 'CM_Data', 'Emass_inv_CBM_FT.dat')
      
   PATH_VBANDS = os.path.join(ROOT_DIR, 'CM_Data', 'valence_bands.pdf')
   PATH_CBANDS = os.path.join(ROOT_DIR, 'CM_Data', 'conduction_bands.pdf')

   VBMk = np.loadtxt(PATH_VBM, dtype = 'cfloat')
   EMASS_VBMk = np.loadtxt(PATH_EMASS_VBM, dtype = 'cfloat')
   CBMk = np.loadtxt(PATH_CBM, dtype = 'cfloat')
   EMASS_CBMk = np.loadtxt(PATH_EMASS_CBM, dtype = 'cfloat')
   
   kpath = generate_kpath(hspoints, Nk_points)
   gpot = np.real(VBMk[:,:2])
   VBpot = VBMk[:,2]
   emass_field_VB = EMASS_VBMk[:,2]
   CBpot = CBMk[:,2]
   emass_field_CB = EMASS_CBMk[:,2]
   
   Ggrid, hdim = generate_G_grid(gdim)
   V_bands = np.zeros((hdim, len(kpath)))
   C_bands = np.zeros((hdim, len(kpath)))     

   # Bandstructure calculation       
   for i,kpoint in enumerate(kpath):
       
       #valence bands
       Hk = generate_Hk(kpoint, Ggrid, hdim, gpot, VBpot, emass_field_VB, gdim, include_emass_correction)
       eigvals, eigens = np.linalg.eigh(Hk)
       
       for j in range(hdim):
           V_bands[j,i] = np.real(eigvals[j])
        
       #conduction bands    
       Hk = generate_Hk(kpoint, Ggrid, hdim, gpot, CBpot, emass_field_CB, gdim, include_emass_correction)
       eigvals, eigens = np.linalg.eigh(Hk)

       for j in range(hdim):
           C_bands[j,i] = np.real(eigvals[j])   

   # Plotting bands        
   kpoints = range(len(kpath))

   # Valence bands
   for j in range(nbands):
       plt.plot(kpoints, V_bands[hdim-j-1,:],c='black')   
       
   plt.ylabel('Energy (meV)', fontsize = 14, fontfamily = 'sans-serif')    
   plt.xticks(np.arange(0, len(kpath), Nk_points), hsnames, fontsize = 20, fontfamily = 'sans-serif')
   plt.savefig(PATH_VBANDS)
   plt.close()
  
   # Conduction bands
   for j in range(nbands):
       plt.plot(kpoints, C_bands[j,:],c='black')   
       
   plt.ylabel('Energy (meV)', fontsize = 14, fontfamily = 'sans-serif')    
   plt.xticks(np.arange(0, len(kpath), Nk_points), hsnames, fontsize = 20, fontfamily = 'sans-serif')
   plt.savefig(PATH_CBANDS)
   plt.close()

   # plotting chrage distribution of some selected eigenvectors
   if plot_eigen:
       print('plotting charge distribution')
       if os.path.exists(ROOT_DIR + '/CM_Data/Eigenvectors/'):
           shutil.rmtree(ROOT_DIR + '/CM_Data/Eigenvectors/')
           
       os.makedirs(ROOT_DIR + '/CM_Data/Eigenvectors/')
       keys = eigen_kpoints.keys()
       #C_eigens = np.zeros((hdim, hdim, len(keys))) * 1j 
       ftdim = 50    
       for i, key in enumerate(keys):
           
           os.makedirs(ROOT_DIR + '/CM_Data/Eigenvectors/' + key)
           kpoint = eigen_kpoints[key][0] * G1m +  eigen_kpoints[key][1] * G2m
       
           #valence bands
           Hk = generate_Hk(kpoint, Ggrid, hdim, gpot, VBpot, emass_field_VB, gdim, include_emass_correction)
           eigvals, eigens = np.linalg.eigh(Hk)

           for b in range(nbands):
               eigen = eigens[:, hdim - b - 1]              
               ldos = np.zeros((ftdim*3,ftdim*3))
      
               for ii in range(-1,2):
                   for jj in range(-1,2):
  
                       RR = R1m * ii + R2m * jj 
                       for i,j in product(range(ftdim), range(ftdim)):
                           rij = R1m * i / ftdim + R2m * j / ftdim + RR
                 
                           FTrans = 0j
                           for kk, Gvec in enumerate(Ggrid):
                               FTrans += eigen[kk] * np.exp(1j*np.dot(Gvec, rij))
          
                           FTrans = FTrans*np.exp(1j*np.dot(kpoint, rij))

                           ldos[i + (ii + 1)*ftdim, j + (jj + 1)*ftdim] = np.absolute(FTrans)**2 / hdim   

               name = os.path.join(ROOT_DIR, 'CM_Data', 'Eigenvectors', key, 'charge_VB_{}.pdf'.format(b))
               plot_moire(ldos, name,  'moire', dtype = 'ldos', colormap = plt.cm.get_cmap('viridis'),repeat = 1)
 
           #conduction bands
           Hk = generate_Hk(kpoint, Ggrid, hdim, gpot, CBpot, emass_field_CB, gdim, include_emass_correction)
           eigvals, eigens = np.linalg.eigh(Hk)

           for b in range(nbands):
               eigen = eigens[:, b]              
               ldos = np.zeros((ftdim*3,ftdim*3))
      
               for ii in range(-1,2):
                   for jj in range(-1,2):
  
                       RR = R1m * ii + R2m * jj 
                       for i,j in product(range(ftdim), range(ftdim)):
                           rij = R1m * i / ftdim + R2m * j / ftdim + RR
                 
                           FTrans = 0j
                           for kk, Gvec in enumerate(Ggrid):
                               FTrans += eigen[kk] * np.exp(1j*np.dot(Gvec, rij))
          
                           FTrans = FTrans*np.exp(1j*np.dot(kpoint, rij))

                           ldos[i + (ii + 1)*ftdim, j + (jj + 1)*ftdim] = np.absolute(FTrans)**2 / hdim   

               name = os.path.join(ROOT_DIR, 'CM_Data', 'Eigenvectors', key, 'charge_CB_{}.pdf'.format(b))
               plot_moire(ldos, name,  'moire', dtype = 'ldos', colormap = plt.cm.get_cmap('viridis'),repeat = 1)
          
          
   # plotting local density of states
   if plot_ldos:
       print('plotting local density of states')
       if os.path.exists(ROOT_DIR + '/CM_Data/Local_dos/'):
           shutil.rmtree(ROOT_DIR + '/CM_Data/Local_dos/')
           
       os.makedirs(ROOT_DIR + '/CM_Data/Local_dos/')
       
       ftdim = 30  
       kmesh = [i*G1m/nkmesh + j*G2m/nkmesh for i,j in product(range(nkmesh), range(nkmesh))]
       
       ldos = np.zeros((ftdim*3,ftdim*3, nbands))
       for i, kpoint in enumerate(kmesh):

           #valence bands
           Hk = generate_Hk(kpoint, Ggrid, hdim, gpot, VBpot, emass_field_VB, gdim, include_emass_correction)
           eigvals, eigens = np.linalg.eigh(Hk)

           for b in range(nbands):
               eigen = eigens[:, hdim - b - 1]              
      
               for ii in range(-1,2):
                   for jj in range(-1,2):
  
                       RR = R1m * ii + R2m * jj 
                       for i,j in product(range(ftdim), range(ftdim)):
                           rij = R1m * i / ftdim + R2m * j / ftdim + RR
                 
                           FTrans = 0j
                           for kk, Gvec in enumerate(Ggrid):
                               FTrans += eigen[kk] * np.exp(1j*np.dot(Gvec, rij))
          
                           FTrans = FTrans*np.exp(1j*np.dot(kpoint, rij))

                           ldos[i + (ii + 1)*ftdim, j + (jj + 1)*ftdim, b] += np.absolute(FTrans)**2 / hdim   
       
       for b in range(nbands):
           name = os.path.join(ROOT_DIR, 'CM_Data', 'Local_dos', 'total_ldos_VB_{}.pdf'.format(b))
           plot_moire(ldos[:,:,b], name,  'moire', dtype = 'ldos', colormap = plt.cm.get_cmap('viridis'),repeat = 1)
 
           #conduction bands
           Hk = generate_Hk(kpoint, Ggrid, hdim, gpot, CBpot, emass_field_CB, gdim, include_emass_correction)
           eigvals, eigens = np.linalg.eigh(Hk)

           for b in range(nbands):
               eigen = eigens[:, b]              
      
               for ii in range(-1,2):
                   for jj in range(-1,2):
  
                       RR = R1m * ii + R2m * jj 
                       for i,j in product(range(ftdim), range(ftdim)):
                           rij = R1m * i / ftdim + R2m * j / ftdim + RR
                 
                           FTrans = 0j
                           for kk, Gvec in enumerate(Ggrid):
                               FTrans += eigen[kk] * np.exp(1j*np.dot(Gvec, rij))
          
                           FTrans = FTrans*np.exp(1j*np.dot(kpoint, rij))

                           ldos[i + (ii + 1)*ftdim, j + (jj + 1)*ftdim, b] += np.absolute(FTrans)**2 / hdim   
       
       for b in range(nbands):
           name = os.path.join(ROOT_DIR, 'CM_Data', 'Local_dos', 'total_ldos_CB_{}.pdf'.format(b))
           plot_moire(ldos[:,:,b], name,  'moire', dtype = 'ldos', colormap = plt.cm.get_cmap('viridis'),repeat = 1)
 
       
