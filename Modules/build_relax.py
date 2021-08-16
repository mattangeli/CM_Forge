import numpy as np
import shutil
import sys 
from params import *

def create_monolayer_relax_inputs():    
    print(' \n Building the Monolayer relaxation folders \n')
    # read atomic positions and species from relax_input_file
    PATH_RELAX_L1 = os.path.join(ROOT_DIR, input_folder, relax_input_file_L1)
    try:
        with open(PATH_RELAX_L1, 'r') as infile:
             lines_l1 = infile.readlines()
                             
    except FileNotFoundError:
        print('\n Relax input file  relative to layer 1 NOT FOUND \n')
    
    PATH_RELAX_L2 = os.path.join(ROOT_DIR, input_folder, relax_input_file_L2)
    try:
        with open(PATH_RELAX_L2, 'r') as infile:
             lines_l2 = infile.readlines()
                             
    except FileNotFoundError:
        print('\n Relax input file  relative to layer 2 NOT FOUND \n')
    
    
    for i in range(nstrain):
        frac_strain = 1 + (i-int(nstrain/2)) * dstrain
        directory = str((i-int(nstrain/2))*3)
           
        lattK=lattice_matrix * np.array([[frac_strain,frac_strain,1] for i in range(3)])
        lattG=lattice_matrix * np.array([[frac_strain, 1 / frac_strain,1] for i in range(3)])   
              
        folder = os.path.join(ROOT_DIR, 'Monolayer', 'Layer_1', 'K' , directory)
    
        if not os.path.exists(folder):
           os.mkdir(folder)
        PATH_sbatch_input = os.path.join(input_folder, sbatch_file)   
        shutil.copy2(PATH_sbatch_input, folder)
           
        os.chdir(folder)
        with open('relax.in', 'w+') as outfile:
                
            towrite = False
            lines_input = lines_l1
            for l,line in enumerate(lines_input):
                strings = line.split() 
    
                if '&IONS' in strings:
                   towrite = True
    
                if '/' in strings and towrite:   
                   line += '\n'
                   line += ' CELL_PARAMETERS {angstrom}' + '\n'         
                   line += str(lattK[0,0]) + ' ' + str(lattK[0,1]) + ' ' + str(lattK[0,2]) 
                   line += '\n'     
                   line += str(lattK[1,0]) + ' ' + str(lattK[1,1]) + ' ' + str(lattK[1,2])   
                   line += '\n'     
                   line += str(lattK[2,0]) + ' ' + str(lattK[2,1]) + ' ' + str(lattK[2,2])    
                   line += '\n'     
                   towrite = False
    
                outfile.write(line)
        
            os.system('sbatch ' + sbatch_file + '\n')
            os.chdir(ROOT_DIR)
          
        folder = os.path.join(ROOT_DIR, 'Monolayer', 'Layer_2', 'K' , directory)
    
        if not os.path.exists(folder):
           os.mkdir(folder)
           PATH_sbatch_input = os.path.join(input_folder, sbatch_file)   
           shutil.copy2(PATH_sbatch_input, folder)
           
        os.chdir(folder)
        with open('relax.in', 'w+') as outfile:
                
            towrite = False
            lines_input = lines_l2
            for l,line in enumerate(lines_input):
                strings = line.split() 
    
                if '&IONS' in strings:
                   towrite = True
    
                if '/' in strings and towrite:   
                   line += '\n'
                   line += ' CELL_PARAMETERS {angstrom}' + '\n'         
                   line += str(lattK[0,0]) + ' ' + str(lattK[0,1]) + ' ' + str(lattK[0,2]) 
                   line += '\n'     
                   line += str(lattK[1,0]) + ' ' + str(lattK[1,1]) + ' ' + str(lattK[1,2])   
                   line += '\n'     
                   line += str(lattK[2,0]) + ' ' + str(lattK[2,1]) + ' ' + str(lattK[2,2])    
                   line += '\n'     
                   towrite = False
    
                outfile.write(line)
        
            os.system('sbatch ' + sbatch_file + '\n')
            os.chdir(ROOT_DIR)
    
        folder = os.path.join(ROOT_DIR, 'Monolayer', 'Layer_1', 'G' , directory)
    
        if not os.path.exists(folder):
           os.mkdir(folder)
           PATH_sbatch_input = os.path.join(input_folder, sbatch_file)   
           shutil.copy2(PATH_sbatch_input, folder)
           
        os.chdir(folder)
        with open('relax.in', 'w+') as outfile:
                
            towrite = False
            lines_input = lines_l1
            for l,line in enumerate(lines_input):
                strings = line.split() 
    
                if '&IONS' in strings:
                   towrite = True
    
                if '/' in strings and towrite:   
                   line += '\n'
                   line += ' CELL_PARAMETERS {angstrom}' + '\n'         
                   line += str(lattG[0,0]) + ' ' + str(lattG[0,1]) + ' ' + str(lattG[0,2]) 
                   line += '\n'     
                   line += str(lattG[1,0]) + ' ' + str(lattG[1,1]) + ' ' + str(lattG[1,2])   
                   line += '\n'     
                   line += str(lattG[2,0]) + ' ' + str(lattG[2,1]) + ' ' + str(lattG[2,2])    
                   line += '\n'     
                   towrite = False
    
    
                outfile.write(line)
        
            os.system('sbatch ' + sbatch_file + '\n')
            os.chdir(ROOT_DIR)
    
        folder = os.path.join(ROOT_DIR, 'Monolayer', 'Layer_2', 'G' , directory)
    
        if not os.path.exists(folder):
           os.mkdir(folder)
           PATH_sbatch_input = os.path.join(input_folder, sbatch_file)   
           shutil.copy2(PATH_sbatch_input, folder)
           
        os.chdir(folder)
        with open('relax.in', 'w+') as outfile:
                
            towrite = False
            lines_input = lines_l2
            for l,line in enumerate(lines_input):
                strings = line.split() 
    
                if '&IONS' in strings:
                   towrite = True
    
                if '/' in strings and towrite:   
                   line += '\n'
                   line += ' CELL_PARAMETERS {angstrom}' + '\n'         
                   line += str(lattG[0,0]) + ' ' + str(lattG[0,1]) + ' ' + str(lattG[0,2]) 
                   line += '\n'     
                   line += str(lattG[1,0]) + ' ' + str(lattG[1,1]) + ' ' + str(lattG[1,2])   
                   line += '\n'     
                   line += str(lattG[2,0]) + ' ' + str(lattG[2,1]) + ' ' + str(lattG[2,2])    
                   line += '\n'     
                   towrite = False
    
                outfile.write(line)
        
            os.system('sbatch ' + sbatch_file + '\n')
            os.chdir(ROOT_DIR)
          
def create_bilayer_relax_inputs():          
    print(' \n Building the Bilayer relaxation folders \n')
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
        print('\n Bilayer relax input file NOT FOUND \n')
    
    
    for i in range(ngrid):
        for j in range(ngrid):
    
            directory = str(i) + '_' + str(j)
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
    
                           if k > 3:
                              atom[1:] = pos[k-1,:] + i*np.asarray([1./ngrid,0.,0.]) + j*np.asarray([0.,1./ngrid,0.])
                           else:
                              atom[1:] = pos[k-1,:] 
    
                           lines_input[l+k] = ''
                           for s in atom:
                               lines_input[l+k] += str(s) + ' '
                           
                           lines_input[l+k] += ' 0 0 1 ' + '\n' 
    #
                     outfile.write(line)
    #             
            os.system('sbatch ' + sbatch_file + '\n')
            os.chdir(ROOT_DIR)
    
    
    
    
    
    
    
