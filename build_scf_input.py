import numpy as np
import shutil
import sys 
from params import *

ROOT_DIR = checkfolders()

# reading some parameters from scf_input_file
try:
    PATH_SCF = os.path.join(ROOT_DIR, input_folder, scf_input_file)
    PATH_NSCF = os.path.join(ROOT_DIR, input_folder, nscf_input_file)    
    PATH_SBATCH = os.path.join(ROOT_DIR, input_folder, sbatch_file)    
    PATH_BANDSX = os.path.join(ROOT_DIR, input_folder, bandsx_file)    
        
    with open(PATH_SCF, 'r') as infile:
         lines_scf = infile.readlines()
    with open(PATH_NSCF, 'r') as infile:
         lines_nscf = infile.readlines()
    with open(PATH_SBATCH, 'r') as infile:
         lines_sbatch = infile.readlines()  
    with open(PATH_BANDSX, 'r') as infile:
         lines_bandsx = infile.readlines()       

except FileNotFoundError:
    sys.exit('SCF, NSCF, SBATCH or BANDSX input file missing in ' + input_folder)         

# reading the relaxed atomic positions and generating the scf input file
for i in range(ngrid):
    for j in range(ngrid):

        directory = str(i) + str(j)
        PATH_Bilayer = os.path.join(ROOT_DIR, 'Bilayer', directory)
        os.chdir(PATH_Bilayer)

# reading relaxed structure
        with open('relax.out', 'r') as infile:
             lines_relax = infile.readlines()
               
             pos = np.zeros([n_atoms, 3])
             species = []
             
             for l,line in enumerate(lines_relax):
                 if line == 'convergence NOT achieved after 100 iterations: stopping':
                      sys.exit('convergence NOT achieved during relax in folder ' + directory)         
                 strings = line.split()
              
                 if 'Begin' in strings:
                     for k in range(3, n_atoms + 3):
                         atom = lines_relax[l+k]
                         species.append(atom.split()[0])
                         pos[k-3,:] = atom.split()[1:]
                 
                     break

# generate scf  and nscf input             
        with open('scf.in', 'w+') as outfile:
             
            lines_input = lines_scf[:]
            for l,line in enumerate(lines_input):
                 strings = line.split()

                 if 'ATOMIC_POSITIONS' in strings:
                     for k in range(1, n_atoms + 1):
                         atom = lines_input[l+k].split()
                         atom[0] = species[k-1]
                         atom[1:] = pos[k-1,:] 

                         lines_input[l+k] = ''
                         for s in atom:
                             lines_input[l+k] += str(s) + ' '
                         lines_input[l+k] += '\n' 

                 outfile.write(line)

        # nscf input 
        with open('nscf.in', 'w+') as outfile:
             
            lines_input = lines_nscf[:]
            for l,line in enumerate(lines_input):
                 strings = line.split()

                 if 'ATOMIC_POSITIONS' in strings:
                     for k in range(1, n_atoms + 1):
                         atom = lines_input[l+k].split()
                         atom[0] = species[k-1]
                         atom[1:] = pos[k-1,:] 

                         lines_input[l+k] = ''
                         for s in atom:
                             lines_input[l+k] += str(s) + ' '
                         lines_input[l+k] += '\n' 
               
                 if 'K_POINTS' in strings:
                    line = 'K_POINTS {crystal} \n'
                    lines_input[l+1] = ' 1  \n'
                    
                    lines_input[l+2] = ''
                    for k in range(len(K_edge)):
                        lines_input[l+2] += str(K_edge[k]) + ' ' 
                    lines_input[l+2] +=  ' 1.0 \n'   
                       
                 outfile.write(line)

        # emass input   
        with open('bands_emass1.in', 'w+') as outfile:

            # effective mass 
            lines_input = lines_nscf[:]
            for l,line in enumerate(lines_input):
                 strings = line.split()
               
                 if '&CONTROL' in strings:
                    lines_input[l+1] = 'calculation = \'bands\' \n'

                 if 'ATOMIC_POSITIONS' in strings:
                     for k in range(1, n_atoms + 1):
                         atom = lines_input[l+k].split()
                         atom[0] = species[k-1]
                         atom[1:] = pos[k-1,:] 

                         lines_input[l+k] = ''
                         for s in atom:
                             lines_input[l+k] += str(s) + ' '
                         lines_input[l+k] += '\n' 
                 
                 gvec = np.asarray([1.0, 0.0, 0.0])        
                 if 'K_POINTS' in strings:
                 
                    line = 'K_POINTS {tpiba_b} \n'
                    lines_input[l+1] = ' 2 \n '
                   
                    lines_input[l+2]= ''
                    for k in range(len(K_edge)):
                        lines_input[l+2] += str(K_edge[k]) + ' '
                    lines_input[l+2] += str(nk_emass) + '\n'
                    
                    lines_input[l+3] = ''
                    for k in range(len(K_edge)):
                        lines_input[l+3] += str(K_edge[k] + emass_frac * gvec[k]) + ' '
                    lines_input[l+3] += str(nk_emass) + '\n'

                 outfile.write(line)


        with open('bands_emass2.in', 'w+') as outfile:
             
            # effective mass 
            lines_input = lines_nscf[:]
            for l,line in enumerate(lines_input):
                 strings = line.split()
                 
                 if '&CONTROL' in strings:
                    lines_input[l+1] = 'calculation = \'bands\' \n'
                 
                 if 'ATOMIC_POSITIONS' in strings:
                     for k in range(1, n_atoms + 1):
                         atom = lines_input[l+k].split()
                         atom[0] = species[k-1]
                         atom[1:] = pos[k-1,:] 

                         lines_input[l+k] = ''
                         for s in atom:
                             lines_input[l+k] += str(s) + ' '
                         lines_input[l+k] += '\n' 
                 
                 gvec = np.asarray([0.0, 1.0, 0.0])        
                 if 'K_POINTS' in strings:
                 
                    line = 'K_POINTS {tpiba_b} \n'
                    lines_input[l+1] = ' 2 \n '
                   
                    lines_input[l+2]= ''
                    for k in range(len(K_edge)):
                        lines_input[l+2] += str(K_edge[k]) + ' '
                    lines_input[l+2] += str(nk_emass) + '\n'
                    
                    lines_input[l+3] = ''
                    for k in range(len(K_edge)):
                        lines_input[l+3] += str(K_edge[k] + emass_frac * gvec[k]) + ' '
                    lines_input[l+3] += str(nk_emass) + '\n'

                 outfile.write(line)

        with open('bandsx1.in', 'w+') as outfile:
        
            lines_input = lines_bandsx[:]    
            for l,line in enumerate(lines_input):
                 strings = line.split()
                 
                 if 'filband' in strings:
                     line = 'filband = \'emass_band1\'' + '\n'
                      
                 outfile.write(line)
                      
        with open('bandsx2.in', 'w+') as outfile:
        
            lines_input = lines_bandsx[:]    
            for l,line in enumerate(lines_input):
                 strings = line.split()
                 
                 if 'filband' in strings:
                     line = 'filband = \'emass_band2\'' + '\n'
                      
                 outfile.write(line)
                                      
        # sbatch file                 
        with open(sbatch_file, 'w+') as outfile:
        
            lines_input = lines_sbatch[:]    
            for l,line in enumerate(lines_input):
                 strings = line.split()
                 
                 if 'srun' in strings:
                     line = 'srun pw.x -input scf.in > scf.out' + '\n'
                     line += 'srun pw.x -input nscf.in > nscf_edge.out ' + '\n'    
                     line += 'srun pw.x -input bands_emass1.in > bands_emass1.out' + '\n' 
                     line += 'srun bands.x -input bandsx1.in > bandsx_emass1.out' + '\n' 
                     line += 'srun pw.x -input bands_emass2.in > bands_emass2.out' + '\n' 
                     line += 'srun bands.x -input bandsx2.in > bandsx_emass2.out' + '\n'                       
                 outfile.write(line)
                          
        os.system('sbatch ' + sbatch_file + '\n')
        os.chdir(ROOT_DIR)






