################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import os
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
from scipy.signal import lfilter, filtfilt

################################################################################
# ZACROS OUTPUTS
################################################################################

class ZacrosOutputs():

    def __init__(self, dirname = '.'):

        self.dirname = dirname

    # -------------------------------------------------------------------
    #  READ GENERAL OUTPUT
    # -------------------------------------------------------------------

    def read_general_output(self, filename = None):

        if filename is None:
            filename = os.path.join(self.dirname, 'general_output.txt')

        keys_vect = ['Temperature'                ,
                     'Pressure'                   ,
                     'Number of gas species'      ,
                     'Gas species names'          ,
                     'Gas species energies'       ,
                     'Gas species molar fractions',
                     'Number of surface species'  ,
                     'Surface species names'      ,
                     'Surface species dentation'  ,
                     'Number of lattice sites'    ,
                     'Number of site types'       ,
                     'Number of clusters'         ,
                     'Number of elementary step'  ]

        fileobj = open(filename, 'rU')
        lines = fileobj.readlines()

        var_vect = []
        for line in lines:
            for key in [ key for key in keys_vect if key in line ]:
                var = line.split()[len(key.split()):]
                for i in range(len(var)):
                    if '.' in var[i]:
                        try: var[i] = float(var[i])
                        except: pass
                    else:
                        try: var[i] = int(var[i])
                        except: pass
                if len(var) == 1:
                    var = var[0]
                var_vect.append(var)

        self.temperature          = var_vect[0]
        self.pressure             = var_vect[1]
        self.n_gas_species        = var_vect[2]
        self.gas_species_names    = var_vect[3]
        self.gas_species_energies = var_vect[4]
        self.gas_mol_fracts       = var_vect[5]
        self.n_surf_species       = var_vect[6]
        self.surf_species_names   = var_vect[7]
        self.surf_species_dent    = var_vect[8]
        self.n_lattice_sites      = var_vect[9]
        self.n_site_types         = var_vect[10]
        self.n_clusters           = var_vect[11]
        self.n_reactions          = var_vect[12]

        self.pre_exponentials = np.zeros(self.n_reactions)

        for n, line in enumerate(lines):
            if 'Reaction network' in line:
                num = n
        
        for n, line in enumerate(lines[num+2:num+2+self.n_reactions]):
            self.pre_exponentials[n] = line.split()[4]

        fileobj.close()

    # -------------------------------------------------------------------
    #  READ SPECNUM OUTPUT
    # -------------------------------------------------------------------

    def read_specnum_output(self, filename = None):

        if filename is None:
            filename = os.path.join(self.dirname, 'specnum_output.txt')

        fileobj = open(filename, 'rU')
        lines = fileobj.readlines()

        n_config = len(lines)-1
        variables = lines[0].split()

        surf_species_names = [ var for var in variables if "*" in var ]
        n_surf_species = len(surf_species_names)

        gas_species_names = variables[5+n_surf_species:]
        n_gas_species = len(gas_species_names)

        n_events       = np.zeros(n_config, dtype = int)
        time           = np.zeros(n_config, dtype = float)
        temperature    = np.zeros(n_config, dtype = float)
        energy         = np.zeros(n_config, dtype = float)
        surf_species   = np.zeros((n_config, n_surf_species), dtype = int)
        surf_coverages = np.zeros((n_config, n_surf_species), dtype = float)
        gas_species    = np.zeros((n_config, n_gas_species), dtype = int)
        gas_prod_rates = np.zeros((n_config, n_gas_species), dtype = float)

        for n, line in enumerate(lines[1:]):
            entries = line.split()
            n_events[n], time[n], temperature[n], energy[n] = entries[1:5]
            for i in range(n_surf_species):
                surf_species[n][i] = entries[5+i]
                surf_coverages[n][i] = surf_species[n][i]/self.n_lattice_sites
            for i in range(n_gas_species):
                gas_species[n][i] = entries[5+n_surf_species+i]

        for n in range(1, n_config):
            dt = time[n] - time[n-1]
            for i in range(n_gas_species):
                gas_prod_rates[n][i] = (gas_species[n][i] - 
                                        gas_species[n-1][i])/dt

        self.surf_species_names  = surf_species_names
        self.gas_species_names   = gas_species_names
        self.n_events_specnum    = n_events
        self.time_specnum        = time
        self.temperature_specnum = temperature
        self.energy_specnum      = energy
        self.surf_species        = surf_species
        self.surf_coverages      = surf_coverages
        self.gas_species         = gas_species
        self.gas_prod_rates      = gas_prod_rates

        fileobj.close()

    # -------------------------------------------------------------------
    #  READ LATTICE OUTPUT
    # -------------------------------------------------------------------

    def read_lattice_output(self, filename = None):

        if filename is None:
            filename = os.path.join(self.dirname, 'lattice_output.txt')

        fileobj = open(filename, 'rU')
        lines = fileobj.readlines()

        n_sites_tot = len(lines)-2
        n_neigh_max = len(lines[0].split())-5

        positions        = np.zeros((n_sites_tot, 2), dtype = float)
        site_types       = np.zeros(n_sites_tot, dtype = int)
        n_site_neighbors = np.zeros(n_sites_tot, dtype = int)
        site_neighbors   = np.empty((n_sites_tot, n_neigh_max), dtype = int)

        for n, line in enumerate(lines[2:]):
            entries = line.split()
            for p in range(2):
                positions[n][p] = entries[1+p]
            site_types[n], n_site_neighbors[n] = entries[3:5]
            for i in range(n_site_neighbors[n]):
                site_neighbors[n][i] = entries[5+i]

        self.positions        = positions
        self.site_types       = site_types
        self.n_site_neighbors = n_site_neighbors
        self.site_neighbors   = site_neighbors

        fileobj.close()

    # -------------------------------------------------------------------
    #  READ PROCSTAT OUTPUT
    # -------------------------------------------------------------------

    def read_procstat_output(self, filename = None):

        if filename is None:
            filename = os.path.join(self.dirname, 'procstat_output.txt')

        fileobj = open(filename, 'rU')
        lines = fileobj.readlines()
        configs_lines = [ line for line in lines if 'configuration' in line ]
        n_configs = len(configs_lines)

        reaction_names = lines[0].split()
        n_reactions = len(reaction_names)

        steps = np.zeros(n_configs, dtype = int)
        time  = np.zeros(n_configs)
        reaction_rates = np.zeros((n_configs, n_reactions), dtype = float)
        reaction_steps = np.zeros((n_configs, n_reactions), dtype = int)

        m = 0
        for n, line in enumerate(lines[1:]):
            if n % 3 == 0:
                steps[m], time[m] = line.split()[2:4]
            elif (n-1) % 3 == 0:
                reaction_rates[m][:] = line.split()[:]
            else:
                reaction_steps[m][:] = line.split()[:]
                m += 1

        self.reaction_names = reaction_names
        self.steps_procstat = steps
        self.time_procstat  = time
        self.reaction_rates = reaction_rates
        self.reaction_steps = reaction_steps

        fileobj.close()

    # -------------------------------------------------------------------
    #  READ HYSTORY OUTPUT
    # -------------------------------------------------------------------

    def read_history_output(self, filename = None):

        if filename is None:
            filename = os.path.join(self.dirname, 'history_output.txt')

        fileobj = open(filename, 'rU')
        lines = fileobj.readlines()

        gas_species_names = lines[0].split()[1:]
        n_gas_species = len(gas_species_names)

        surf_species_names = lines[1].split()[1:]
        surf_species_dict = OrderedDict()
        surf_species_dict[0] = '*'
        for n, spec in enumerate(surf_species_names):
            surf_species_dict[n+1] = spec

        for n, line in enumerate(lines):
            if 'configuration' in line:
                start = n
                break

        for n, line in enumerate(lines[start+1:]):
            if 'configuration' in line:
                n_sites = n-1
                break

        configs_lines = [ line for line in lines if 'configuration' in line ]
        n_configs = len(configs_lines)

        steps       = np.zeros(n_configs, dtype = int)
        time        = np.zeros(n_configs, dtype = float)
        temperature = np.zeros(n_configs, dtype = float)
        energy      = np.zeros(n_configs, dtype = float)
        adsorbate_number    = np.zeros((n_configs, n_sites), dtype = int)
        surf_species_config = np.zeros((n_configs, n_sites), dtype = int)
        dentation_config    = np.zeros((n_configs, n_sites), dtype = int)
        gas_species_prod    = np.zeros((n_configs, n_gas_species), dtype = int)

        for n, line in enumerate(configs_lines):
            steps[n], time[n], temperature[n], energy[n] = line.split()[2:6]

        for n in range(n_configs):
            for s in range(n_sites):
                entries = lines[start+1+n*(n_sites+2)+s].split()
                adsorbate_number[n][s]    = entries[1]
                surf_species_config[n][s] = entries[2]
                dentation_config[n][s]    = entries[3]
            production_vect = lines[start+1+n*(n_sites+2)+n_sites].split()
            for g in range(n_gas_species):
                gas_species_prod[n][g] = int(production_vect[g])

        self.surf_species_dict   = surf_species_dict
        self.steps_history       = steps
        self.time_history        = time
        self.temperature_history = temperature
        self.energy_history      = energy
        self.adsorbate_number    = adsorbate_number
        self.surf_species_config = surf_species_config
        self.dentation_config    = dentation_config
        self.gas_species_prod    = gas_species_prod

        fileobj.close()

    # -------------------------------------------------------------------
    #  PLOT GAS RATES
    # -------------------------------------------------------------------

    def plot_gas_prod_rates(self, filters = None):

        plot_variables(x       = self.time_specnum               ,
                       y       = self.gas_prod_rates             ,
                       labels  = self.gas_species_names          ,
                       title   = 'gas species production rates'  ,
                       xlabel  = 'time [s]'                      ,
                       ylabel  = 'production rates [molecules/s]',
                       filters = filters                         )

    # -------------------------------------------------------------------
    #  PLOT SURF SPECIES
    # -------------------------------------------------------------------

    def plot_surf_coverages(self, filters = None):

        plot_variables(x       = self.time_specnum          ,
                       y       = self.surf_coverages        ,
                       labels  = self.surf_species_names    ,
                       title   = 'surface species coverages',
                       xlabel  = 'time [s]'                 ,
                       ylabel  = 'coverages [-]'            ,
                       ymin    = 0.                         ,
                       ymax    = 1.                         ,
                       filters = filters                    )

    # -------------------------------------------------------------------
    #  PLOT REACTION RATES
    # -------------------------------------------------------------------

    def plot_reaction_rates(self, filters = None):

        plot_variables(x       = self.time_procstat  ,
                       y       = self.reaction_rates ,
                       labels  = self.reaction_names ,
                       title   = 'reaction rates'    ,
                       xlabel  = 'time [s]'          ,
                       ylabel  = 'reaction rates [-]',
                       ymin    = 0.                  ,
                       filters = filters             )

    # -------------------------------------------------------------------
    #  CHECK ALL CONVERGED
    # -------------------------------------------------------------------

    def check_all_converged(self,
                            delta_gas = 0.1 ,
                            delta_ads = 0.1 ,
                            delta_rxn = 0.1 ,
                            portion   = 0.8 ,
                            pieces    = 4   ,
                            filters   = None):

        convergence_reached = False

        self.read_general_output()
        self.read_specnum_output()
        self.read_procstat_output()

        gas_conv = check_convergence(y       = self.gas_prod_rates,
                                     delta   = delta_gas          ,
                                     portion = portion            ,
                                     pieces  = pieces             ,
                                     filters = filters            )

        ads_conv = check_convergence(y       = self.surf_coverages,
                                     delta   = delta_ads          ,
                                     portion = portion            ,
                                     pieces  = pieces             ,
                                     filters = filters            )

        rxn_conv = check_convergence(y       = self.reaction_rates,
                                     delta   = delta_rxn          ,
                                     portion = portion            ,
                                     pieces  = pieces             ,
                                     filters = filters            )

        if gas_conv is True and ads_conv is True and rxn_conv is True:
            convergence_reached = True

        return convergence_reached

    # -------------------------------------------------------------------
    #  CHECK ALL CONVERGED
    # -------------------------------------------------------------------

    def check_gas_surf_converged(self,
                                 delta_gas = 0.1 ,
                                 delta_ads = 0.1 ,
                                 delta_rxn = 0.1 ,
                                 portion   = 0.8 ,
                                 pieces    = 6   ,
                                 filters   = None):

        convergence_reached = False

        self.read_general_output()
        self.read_specnum_output()
        self.read_procstat_output()

        gas_conv = check_convergence(y       = self.gas_prod_rates,
                                     delta   = delta_gas          ,
                                     portion = portion            ,
                                     pieces  = pieces             ,
                                     filters = filters            )

        ads_conv = check_convergence(y       = self.surf_coverages,
                                     delta   = delta_ads          ,
                                     portion = portion            ,
                                     pieces  = pieces             ,
                                     filters = filters            )

        if gas_conv is True and ads_conv is True:
            convergence_reached = True

        return convergence_reached

    # -------------------------------------------------------------------
    #   GET PRODUCTION RATES
    # -------------------------------------------------------------------
    
    def print_results(self, portion = 0.4, resultsfile = 'RESULTS.txt'):
        
        f = open(resultsfile, 'w+')
        
        self.gas_prod_rates_final = average_results(self.gas_prod_rates,
                                                    portion = portion  )
        self.surf_coverages_final = average_results(self.surf_coverages,
                                                    portion = portion  )
        self.reaction_rates_final = average_results(self.reaction_rates,
                                                    portion = portion  )
        
        col = 34
        
        print('PRE EXPONENTIALS', file = f)
        
        for i in range(len(self.reaction_names)-1):
            print(self.reaction_names[i+1].ljust(col),
                  ' {0:+13.6e}'.format(self.pre_exponentials[i]), file = f)
        
        print('', file = f)
        
        print('REACTION RATES', file = f)
        
        for i in range(len(self.reaction_names)):
            print(self.reaction_names[i].ljust(col),
                  ' {0:+13.6e}'.format(self.reaction_rates_final[i]), file = f)
        
        print('', file = f)
        
        print('GAS SPECIES PRODUCTION RATES', file = f)
        
        for i in range(len(self.gas_species_names)):
            print(self.gas_species_names[i].ljust(col),
                  ' {0:+13.6e}'.format(self.gas_prod_rates_final[i]), file = f)
        
        print('', file = f)
        
        print('SURFACE SPECIES COVERAGES', file = f)
        
        for i in range(len(self.surf_species_names)):
            print(self.surf_species_names[i].ljust(col),
                  ' {0:+13.6e}'.format(self.surf_coverages_final[i]), file = f)
        
        f.close()

################################################################################
# PLOT PLOTS
################################################################################

def plot_variables(x, y, labels, title, xlabel, ylabel, ymin = None,
                   ymax = None, filters = None):

    plt.figure(figsize = (15, 10), dpi = (150))

    plt.rc('font', size = 18)
    plt.rc('axes', titlesize = 18)
    plt.rc('axes', labelsize = 18)
    plt.rc('xtick', labelsize = 14)
    plt.rc('ytick', labelsize = 14)
    plt.rc('legend', fontsize = 10)
    plt.rc('figure', titlesize = 20)

    ax = plt.subplot(111)

    if filters is not None:
        y = reduce_noise(y, filters = filters)

    if ymin is None:
        ymin = np.min(y)*1.1

    if ymax is None:
        ymax = np.max(y)*1.1

    if len(y.shape) == 1:
        ax.plot(x, y, label = labels[0])
    else:
        for i in range(y.shape[1]):
            ax.plot(x, y[:,i], label = labels[i])

    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.axis([min(x), max(x), ymin, ymax])
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))

    plt.savefig(title+'.png')
    plt.close()

################################################################################
# REDUCE NOISE
################################################################################

def reduce_noise(y, filters = 15):

    y_filtered = np.copy(y)
    b = [1./filters]*filters
    a = 1

    if len(y.shape) == 1:
        y_filtered = filtfilt(b, a, y)
    else:
        for i in range(y.shape[1]):
            y_filtered[:,i] = filtfilt(b, a, y[:,i])

    return y_filtered

################################################################################
# CHECK CONVERGENCE
################################################################################

def check_convergence(y, delta = 0.05, portion = 0.8, pieces = 4, filters = 15):

    try:
        if filters is not None:
            y = reduce_noise(y, filters = filters)
    except: pass

    if len(y.shape) == 1:

        y_cut = y[int((1-portion)*len(y)):]

        check = convergence(y_cut = y_cut, delta = delta, pieces = pieces)

    else:

        for i in range(y.shape[1]):

            y_cut = y[:,i][int((1-portion)*len(y[:,i])):]

            check = convergence(y_cut = y_cut, delta = delta, pieces = pieces)

            if check is False:
                break

    return check

################################################################################
# CONVERGENCE
################################################################################

def convergence(y_cut, delta = 0.05, pieces = 4):

    check = False

    if len(y_cut) > pieces:

        y_split = np.array_split(y_cut, pieces)
        y_ave = np.zeros(pieces)
    
        for i in range(len(y_split)):
            y_ave[i] = np.sum(y_split[i])/len(y_split[i])
        
        if np.max(np.abs(y_ave)) > 0.:
            diff = abs(np.max(y_ave)-np.min(y_ave))/np.max(np.abs(y_ave))
        else:
            diff = 0.

        if diff < delta:
            check = True

    return check

################################################################################
# AVERAGE RESULTS
################################################################################

def average_results(y, portion = 0.8):

    y_ave = np.zeros(y.shape[1])
    
    for i in range(y.shape[1]):
        y_cut = y[:,i][int((1-portion)*len(y[:,i])):]
        y_ave[i] = np.sum(y_cut)/len(y_cut)

    return y_ave

################################################################################
# END
################################################################################
