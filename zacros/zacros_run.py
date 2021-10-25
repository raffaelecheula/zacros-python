################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import os, subprocess, time
from zacros_outputs import ZacrosOutputs
from zacros_utils import remove_zacros_outputs
import zacros_globals

################################################################################
# ZACROS GLOBALS INITIALIZE
################################################################################

zacros_globals.initialize()

################################################################################
# ZACROS RUN TO CONVERGENCE
################################################################################

def zacros_run_to_convergence(zacros_exe   = './zacros.x' , 
                              dirname      = '.'          ,
                              from_scratch = True         ,
                              resultsfile  = 'RESULTS.txt', 
                              delta_gas    = 0.1          ,
                              delta_ads    = 0.1          ,
                              delta_rxn    = 0.1          ,
                              portion      = 0.4          ,
                              pieces       = 4            ,
                              plot_plots   = True         ,
                              check_conv   = 'all'        ,
                              filters      = None         ):

    workdir = os.getcwd()
    os.chdir(dirname)

    if from_scratch is True:
        remove_zacros_outputs()

    if check_conv not in ('all', 'gas surf'):
        raise ValueError

    converged = False
    
    while converged is False:
        
        subprocess.call(zacros_exe)
        
        zout = ZacrosOutputs()

        if check_conv == 'all':
            converged = zout.check_all_converged(delta_gas = delta_gas,
                                                 delta_ads = delta_ads,
                                                 delta_rxn = delta_rxn,
                                                 portion   = portion  ,
                                                 pieces    = pieces   ,
                                                 filters   = filters  )
        elif check_conv == 'gas surf':
            converged = zout.check_gas_surf_converged(delta_gas = delta_gas,
                                                      delta_ads = delta_ads,
                                                      portion   = portion  ,
                                                      pieces    = pieces   ,
                                                      filters   = filters  )

    zout.print_results(portion = portion, resultsfile = resultsfile)

    if plot_plots is True:
        zout.plot_gas_prod_rates()
        zout.plot_surf_coverages()
        zout.plot_reaction_rates()

    with zacros_globals.number.get_lock():
        zacros_globals.number.value -= 1

    os.chdir(workdir)

################################################################################
# END
################################################################################
