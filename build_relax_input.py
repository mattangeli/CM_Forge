import numpy as np
import shutil
import sys 
from params import *

ROOT_DIR = checkfolders()

# read atomic positions and species from relax_input_file
try:
    PATH_RELAX = os.path.join(ROOT_DIR, input_folder, relax_input_file)
    with open(PATH_RELAX, 'r') as infile:
         lines = infile.readlines()
         pos = np.zeros([n_atoms, 3])
         species = []
         
         for l,line in enumerate(lines):
             strings = line.split()
          
             if 'ATOMIC_POSITIONS' in strings:
                 for i in range(1, n_atoms + 1):
                     atom = lines[l+i]
                     species.append(atom.split()[0])
                     pos[i-1,:] = atom.split()[1:]
             
                 break
                
except FileNotFoundError:
    print('\n Relax input file NOT FOUND \n')


for i in range(ngrid):
    for j in range(ngrid):

        directory = str(i) + str(j)
        PATH_Bilayer = os.path.join(ROOT_DIR, 'Bilayer', directory)
        if not os.path.exists(PATH_Bilayer):
           os.mkdir(PATH_Bilayer)
        PATH_sbatch_input = os.path.join(input_folder, sbatch_file)   
        shutil.copy2(PATH_sbatch_input, PATH_Bilayer)
        os.chdir(PATH_Bilayer)

        with open('relax.in', 'w+') as outfile:
             
             lines_input = lines
             for l,line in enumerate(lines_input):
                 strings = line.split()
#
                 if 'ATOMIC_POSITIONS' in strings:
                    for k in range(1, n_atoms + 1):
                        atom = lines_input[l+k].split()
                        atom[0] = species[k-1]
                        atom[1:] = pos[k-1,:] + i*np.asarray([1./ngrid,0.,0.]) + j*np.asarray([0.,1./ngrid,0.])
#                        
                        lines_input[l+k] = ''
                        for s in atom:
                            lines_input[l+k] += str(s) + ' '
                        lines_input[l+k] += '\n' 
#
                 outfile.write(line)
#             
        os.system('sbatch ' + sbatch_file + '\n')
        os.chdir(ROOT_DIR)









