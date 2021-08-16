import params
from Modules.build_relax import *
from Modules.build_electronic_structure import *
from Modules.Relaxation import relaxation
from Modules.read_moire_potential import *
from Modules.CM_bands import *

# generate and execute QE DFT calculations
generate_monolayer_relax = False
generate_bilayer_relax = False
generate_bilayer_bands = False

# Analize data, relax and compute CM bands
relax_moire_cell = 1
read_moire_potential = 1
include_relaxation = 1
CM_bands = 1
include_emass_correction = 1
plot_eigen = 1
plot_ldos = 0

if generate_monolayer_relax:
   
   create_monolayer_relax_inputs()

if generate_bilayer_relax:

   create_bilayer_relax_inputs()

if relax_moire_cell:
   
   relaxation()

if generate_bilayer_bands:
   
   generate_bilayer_electronic_structure() 

if read_moire_potential:
   
   read_fields(include_relaxation)

if CM_bands:
   
   continuum_model_bands(include_emass_correction, plot_eigen, plot_ldos)












































