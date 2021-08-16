import numpy as np
import shutil
import sys 
import params as par
from utils import *

ROOT_DIR = checkfolders()

def relaxation():
    print('The moirè lenght is:', np.linalg.norm(par.R1m), 'A')
    # read atomic positions and species from relax_input_file
    PATH_RELAX_L1 = os.path.join(ROOT_DIR, 'Monolayer', 'Layer_1')
    PATH_RELAX_L2 = os.path.join(ROOT_DIR, 'Monolayer', 'Layer_2')
    
    strain_L1, Kt, Gt = strainanalysis(PATH_RELAX_L1)
    strain_L2, Kb, Gb = strainanalysis(PATH_RELAX_L2)        
                    
    PATH_RELAX_GSFE = os.path.join(ROOT_DIR, 'Bilayer')
    GSFE, GSFE_k, GSFEdic = readGSFE(PATH_RELAX_GSFE)
    
    plot_moire(GSFE, ROOT_DIR + '/CM_Data/GSFE_unrelax.pdf', 'unit')
    
    u_rel = relax(GSFE_k, Kt, Gt, Kb, Gb)
    
    GSFE_relax, GSFEk_relax = include_relax(GSFE_k, u_rel)    
    plot_moire(GSFE_relax, ROOT_DIR + '/CM_Data/GSFE_relax.pdf', 'moire', dtype = 'field', colormap = plt.cm.get_cmap('bwr') )
    
    
