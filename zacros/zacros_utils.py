################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import numpy as np
import os, shutil
from itertools import product

################################################################################
# COPY INPUT FILES
################################################################################

def copy_zacros_inputs(dirname = '.'):

    if not os.path.isfile(os.path.join(dirname, 'simulation_input.dat')):
        try: shutil.copy2('simulation_input.dat', dirname)
        except: pass
    if not os.path.isfile(os.path.join(dirname, 'energetics_input.dat')):
        try: shutil.copy2('energetics_input.dat', dirname)
        except: pass
    if not os.path.isfile(os.path.join(dirname, 'mechanism_input.dat')):
        try: shutil.copy2('mechanism_input.dat', dirname)
        except: pass
    if not os.path.isfile(os.path.join(dirname, 'lattice_input.dat')):
        try: shutil.copy2('lattice_input.dat', dirname)
        except: pass

################################################################################
# MOLAR FRACTS VECTORS
################################################################################

def get_molar_fracts_vectors(n_gas_specs, x_min_vect, deltax = 0.01, 
                             epsi = 1e-3):

    v = np.arange(0., 1.+deltax, deltax)

    w = [ t for t in product(v, repeat = n_gas_specs) if sum(t) < 1.+epsi 
          if np.all(np.array(t) > np.array([i-epsi for i in x_min_vect])) ]

    return np.array(w)

################################################################################
# MOLAR FRACTS DIRNAME
################################################################################

def molar_fracts_dirname(gas_species_names, gas_molar_fracs):

    dirname = ''
    for i in range(len(gas_species_names)):
        dirname += '_{0}_{1:3.2f}'.format(gas_species_names[i],
                                         gas_molar_fracs[i])

    return dirname

################################################################################
# REMOVE OUTPUT FILES
################################################################################

def remove_zacros_outputs(dirname = '.'):

    try: os.remove(os.path.join(dirname, 'general_output.txt'))
    except: pass
    try: os.remove(os.path.join(dirname, 'history_output.txt'))
    except: pass
    try: os.remove(os.path.join(dirname, 'lattice_output.txt'))
    except: pass
    try: os.remove(os.path.join(dirname, 'procstat_output.txt'))
    except: pass
    try: os.remove(os.path.join(dirname, 'specnum_output.txt'))
    except: pass
    try: os.remove(os.path.join(dirname, 'restart.inf'))
    except: pass

################################################################################
# END
################################################################################
