#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import os, time
import numpy as np
import multiprocessing, time
from zacros_run import zacros_run_to_convergence
from zacros_inputs import (ZacrosInputs, Cluster, Step, ClusterVariant,
                           StepVariant)
from zacros_utils import (copy_zacros_inputs, get_molar_fracts_vectors,
                          molar_fracts_dirname)
import zacros_globals

################################################################################
# INITIALIZE ZACROS INPUTS
################################################################################

zin = ZacrosInputs()

zin.random_seed        = 154686

zin.temperature        = 100.0 # [K]
zin.ramp               = None # [K/s]
zin.pressure           = 2.0 # [bar]

zin.gas_specs_names    = [    'A',     'B',     'C']
zin.gas_energies       = [  0.000,   0.000,   0.000] # [eV]
zin.gas_molec_weights  = [  1.000,   1.000,   1.000] # [uma]
zin.gas_molar_fracs    = [  0.500,   0.500,   0.000] # [-]

zin.surf_specs_names   = ['A*', 'B*']

zin.acquirement        = 'on event'
zin.snapshots          = 1000 # to history_output.dat
zin.process_statistics = 1000 # to procstat_output.dat
zin.species_numbers    = 1000 # to specnum_output.dat

zin.event_report       = False

zin.max_steps          = None # [s]
zin.max_time           = None # [s]

zin.wall_time          = 60 # [s]

zin.restart            = True

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
parameterfile = 'PARAMETERS.txt'
n_max_process = 8

workdir = os.getcwd()
zacros_globals.initialize()

k1 = 1e0

vec_ads  = [1e0]
vec_diff = [0e0]
vec_span = [1e-4, 1e-2, 1e0, 1e2, 1e4]

sim_num = 0

for (k2, k3, k4, k5, k6, k7) in [ (k2, k3, k4, k5, k6, k7)
                                   for k2 in vec_ads
                                   for k3 in vec_ads
                                   for k4 in vec_ads
                                   for k5 in vec_span
                                   for k6 in vec_diff 
                                   for k7 in vec_diff ]:

    sim_num += 1

    while zacros_globals.number.value >= n_max_process:
        time.sleep(10)

    dirname = 'sumulation_{:04d}'.format(sim_num)

    if os.path.isfile(os.path.join(dirname, resultsfile)):
        break

    # -------------------------------------------------------------------
    #  ENERGETICS INPUT
    # -------------------------------------------------------------------

    zin.clusters = []

    # Binding Energies of A*
    cls = Cluster(configuration = 'A*')
    cls.add_variant(sites = 'top', eng = -0.0)
    #cls.add_variant(sites = 'brg', eng = -0.0)
    zin.clusters.append(cls)
    
    # Binding Energies of B*
    cls = Cluster(configuration = 'B*')
    cls.add_variant(sites = 'top', eng = -0.0)
    #cls.add_variant(sites = 'brg', eng = -0.0)
    zin.clusters.append(cls)

    # -------------------------------------------------------------------
    #  MECHANISM INPUT
    # -------------------------------------------------------------------

    zin.steps = []

    # Adsorption of A
    stp = Step(reaction = 'A + * <=> A*')
    stp.add_variant(sites = 'top', A_fwd = k1, A_rev = k2, Eact = 0.)
    #stp.add_variant(sites = 'brg', A_fwd = 1., A_rev = 1., Eact = 0.)
    zin.steps.append(stp)
    
    # Adsorption of B
    stp = Step(reaction = 'B + * <=> B*')
    stp.add_variant(sites = 'top', A_fwd = k3, A_rev = k4, Eact = 0.)
    #stp.add_variant(sites = 'brg', A_fwd = 1., A_rev = 1., Eact = 0.)
    zin.steps.append(stp)
    
    # Reaction between A* and B*
    stp = Step(reaction = 'A* + B* => * + * + C')
    stp.add_variant(sites = ('top', 'top'), A_fwd = k5, Eact = 0.)
    #stp.add_variant(sites = ('top', 'brg'), A_fwd = 1., Eact = 0.)
    #stp.add_variant(sites = ('brg', 'top'), A_fwd = 1., Eact = 0.)
    #stp.add_variant(sites = ('brg', 'brg'), A_fwd = 1., Eact = 0.)
    zin.steps.append(stp)
    
    # Surface Diffusion of A*
    stp = Step(reaction = 'A* + * => * + A*')
    stp.add_variant(sites = ('top', 'top'), A_fwd = k6, Eact = 0.)
    #stp.add_variant(sites = ('top', 'brg'), A_fwd = 1., Eact = 0.)
    #stp.add_variant(sites = ('brg', 'top'), A_fwd = 1., Eact = 0.)
    #stp.add_variant(sites = ('brg', 'brg'), A_fwd = 1., Eact = 0.)
    zin.steps.append(stp)
    
    # Surface Diffusion of B*
    stp = Step(reaction = 'B* + * => * + B*')
    stp.add_variant(sites = ('top', 'top'), A_fwd = k7, Eact = 0.)
    #stp.add_variant(sites = ('top', 'brg'), A_fwd = 1., Eact = 0.)
    #stp.add_variant(sites = ('brg', 'top'), A_fwd = 1., Eact = 0.)
    #stp.add_variant(sites = ('brg', 'brg'), A_fwd = 1., Eact = 0.)
    zin.steps.append(stp)
    
    try: os.mkdir(dirname)
    except: pass

    copy_zacros_inputs(dirname = dirname)

    os.chdir(dirname)

    zin.print_simulation_input()
    zin.print_energetics_input()
    zin.print_mechanism_input()

    f = open(parameterfile, 'w+')
    for ki in (k1, k2, k3, k4, k5, k6, k7):
        print('{0:+8.4E}'.format(ki), file = f)
    f.close()

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
