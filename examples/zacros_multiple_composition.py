#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import os
import numpy as np
import multiprocessing, time
from zacros_run import zacros_run_to_convergence
from zacros_inputs import ZacrosInputs
from zacros_utils import (copy_zacros_inputs, get_molar_fracts_vectors,
                          molar_fracts_dirname)
import zacros_globals

################################################################################
# INITIALIZE ZACROS INPUTS
################################################################################

zin = ZacrosInputs()

zin.random_seed        = 154686

zin.temperature        = 500.0 # [K]
zin.ramp               = None # [K/s]
zin.pressure           = 1.0 # [bar]

zin.gas_specs_names    = [    'A',     'B',     'C']
zin.gas_energies       = [  0.000,   0.000,   0.000] # [eV]
zin.gas_molec_weights  = [  1.000,   1.000,   1.000] # [uma]
zin.gas_molar_fracs    = [  0.000,   0.000,   0.000] # [-]

zin.surf_specs_names   = ['A*', 'B*']

zin.acquirement        = 'on event'
zin.snapshots          = 1000 # to history_output.dat
zin.process_statistics = 1000 # to procstat_output.dat
zin.species_numbers    = 1000 # to specnum_output.dat

zin.event_report       = False

zin.max_steps          = None # [s]
zin.max_time           = None # [s]

zin.wall_time          = 20 # [s]

zin.restart            = True

zin.print_simulation_input()

################################################################################
# MULTUIPLE ZACROS RUNS
################################################################################

zacros_exe   = zacros_globals.zacros_exe
from_scratch = True
resultsfile  = 'RESULTS.txt'
delta_gas    = 0.2
delta_ads    = 0.1
delta_rxn    = 0.2
pieces       = 6
plot_plots   = True
check_conv   = 'gas surf'
filters      = None

subdirname    = '.'
n_max_process = 8
deltax        = 0.1
x_min_vect    = [0.1, 0.1, 0.0]

workdir = os.getcwd()
zacros_globals.initialize()

n_gas_specs = len(zin.gas_specs_names)
gas_species_names = zin.gas_specs_names

molar_fracts_vectors = get_molar_fracts_vectors(n_gas_specs = n_gas_specs,
                                                deltax      = deltax     ,
                                                x_min_vect  = x_min_vect )

for i in range(len(molar_fracts_vectors)):

    while zacros_globals.number.value >= n_max_process:
        time.sleep(10)

    gas_molar_fracs = molar_fracts_vectors[i]

    dirname = molar_fracts_dirname(gas_species_names, gas_molar_fracs)
    
    try: os.mkdir(dirname)
    except: pass

    if os.path.isfile(os.path.join(dirname, resultsfile)) is False:

        copy_zacros_inputs(dirname = dirname)
    
        os.chdir(dirname)
    
        zin.gas_molar_fracs = gas_molar_fracs
        zin.print_simulation_input()
    
        run = multiprocessing.Process(target = zacros_run_to_convergence,
                                      args   = (zacros_exe  ,
                                                subdirname  ,
                                                from_scratch,
                                                resultsfile ,
                                                delta_gas   ,
                                                delta_ads   ,
                                                delta_rxn   ,
                                                pieces      ,
                                                plot_plots  ,
                                                check_conv  ,
                                                filters     ))
    
        run.start()

        with zacros_globals.number.get_lock():
            zacros_globals.number.value += 1

    os.chdir(workdir)

################################################################################
# END
################################################################################
