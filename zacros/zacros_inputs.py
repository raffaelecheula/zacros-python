################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import os, re, itertools

################################################################################
# ZACROS INPUTS
################################################################################

class ZacrosInputs():

    separation_string = '#' * 80

    def __init__(self,
                 dirname                        = '.',   # str
                 random_seed                    = None,  # int
                 tol_dx_newton                  = None,  # float
                 tol_rhs_newton                 = None,  # float
                 max_newton_iter                = None,  # int
                 n_gauss_pts                    = None,  # int
                 temperature                    = None,  # float
                 ramp                           = None,  # float
                 pressure                       = None,  # float
                 gas_specs_names                = None,  # list(str)
                 gas_energies                   = None,  # list(float)
                 gas_molec_weights              = None,  # list(float)
                 gas_molar_fracs                = None,  # list(float)
                 surf_specs_names               = None,  # list(str)
                 surf_specs_dent                = None,  # list(int)
                 acquirement                    = None,  # str
                 snapshots                      = None,  # tuple | int | float
                 process_statistics             = None,  # tuple | int | float
                 species_numbers                = None,  # tuple | int | float
                 energetics_lists               = None,  # tuple | int | float
                 process_lists                  = None,  # tuple | int | float
                 event_report                   = False, # str | bool
                 max_steps                      = None,  # int
                 max_time                       = None,  # float
                 wall_time                      = None,  # int
                 restart                        = True,  # bool
                 enable_stiffness_scaling       = False, # bool
                 stiffness_scale_all            = False, # bool
                 check_every                    = None,  # int
                 min_separation                 = None,  # float
                 max_separation                 = None,  # float
                 max_qequil_separation          = None,  # float
                 tol_part_equil_ratio           = None,  # float
                 stiffn_coeff_threshold         = None,  # float
                 scaling_factor                 = None,  # float
                 cache_product_clusters         = False, # bool
                 debug_report_global_energetics = False, # bool
                 debug_report_processes         = False, # bool
                 debug_newton_method            = False, # bool
                 debug_check_processes          = False, # bool
                 debug_check_lattice            = False, # bool
                 clusters                       = None,  # list
                 steps                          = None,  # list
                 active_sites                   = None): # list

        if acquirement is None:
            acquirement = ''
        if clusters is None:
            clusters = []
        if steps is None:
            steps = []
        if active_sites is None:
            active_sites = []

        self.dirname                        = dirname
        self.random_seed                    = random_seed                   
        self.tol_dx_newton                  = tol_dx_newton                 
        self.tol_rhs_newton                 = tol_rhs_newton                
        self.max_newton_iter                = max_newton_iter               
        self.n_gauss_pts                    = n_gauss_pts                   
        self.temperature                    = temperature                   
        self.ramp                           = ramp                          
        self.pressure                       = pressure                      
        self.gas_specs_names                = gas_specs_names               
        self.gas_energies                   = gas_energies                  
        self.gas_molec_weights              = gas_molec_weights             
        self.gas_molar_fracs                = gas_molar_fracs               
        self.surf_specs_names               = surf_specs_names              
        self.surf_specs_dent                = surf_specs_dent               
        self.acquirement                    = acquirement                   
        self.snapshots                      = snapshots                     
        self.process_statistics             = process_statistics            
        self.species_numbers                = species_numbers               
        self.energetics_lists               = energetics_lists              
        self.process_lists                  = process_lists                 
        self.event_report                   = event_report                  
        self.max_steps                      = max_steps                     
        self.max_time                       = max_time                      
        self.wall_time                      = wall_time                     
        self.restart                        = restart                       
        self.enable_stiffness_scaling       = enable_stiffness_scaling      
        self.stiffness_scale_all            = stiffness_scale_all           
        self.check_every                    = check_every                   
        self.min_separation                 = min_separation                
        self.max_separation                 = max_separation                
        self.max_qequil_separation          = max_qequil_separation         
        self.tol_part_equil_ratio           = tol_part_equil_ratio          
        self.stiffn_coeff_threshold         = stiffn_coeff_threshold        
        self.scaling_factor                 = scaling_factor                
        self.cache_product_clusters         = cache_product_clusters        
        self.debug_report_global_energetics = debug_report_global_energetics
        self.debug_report_processes         = debug_report_processes        
        self.debug_newton_method            = debug_newton_method           
        self.debug_check_processes          = debug_check_processes         
        self.debug_check_lattice            = debug_check_lattice           
        self.clusters                       = clusters                      
        self.steps                          = steps                         
        self.active_sites                   = active_sites                  

    # -------------------------------------------------------------------
    #  PRINT SIMULATION INPUT
    # -------------------------------------------------------------------

    def print_simulation_input(self, filename = None):

        try: os.mkdir(self.dirname)
        except: pass

        if filename is None:
            filename = os.path.join(self.dirname, 'simulation_input.dat')

        random_seed                    = self.random_seed
        tol_dx_newton                  = self.tol_dx_newton
        tol_rhs_newton                 = self.tol_rhs_newton
        max_newton_iter                = self.max_newton_iter
        n_gauss_pts                    = self.n_gauss_pts
        temperature                    = self.temperature
        ramp                           = self.ramp
        pressure                       = self.pressure
        gas_specs_names                = self.gas_specs_names
        gas_energies                   = self.gas_energies
        gas_molec_weights              = self.gas_molec_weights
        gas_molar_fracs                = self.gas_molar_fracs
        surf_specs_names               = self.surf_specs_names
        surf_specs_dent                = self.surf_specs_dent
        acquirement                    = self.acquirement
        snapshots                      = self.snapshots
        process_statistics             = self.process_statistics
        species_numbers                = self.species_numbers
        energetics_lists               = self.energetics_lists
        process_lists                  = self.process_lists
        event_report                   = self.event_report
        max_steps                      = self.max_steps
        max_time                       = self.max_time
        wall_time                      = self.wall_time
        restart                        = self.restart
        enable_stiffness_scaling       = self.enable_stiffness_scaling
        stiffness_scale_all            = self.stiffness_scale_all
        check_every                    = self.check_every
        min_separation                 = self.min_separation
        max_separation                 = self.max_separation
        max_qequil_separation          = self.max_qequil_separation
        tol_part_equil_ratio           = self.tol_part_equil_ratio
        stiffn_coeff_threshold         = self.stiffn_coeff_threshold
        scaling_factor                 = self.scaling_factor
        cache_product_clusters         = self.cache_product_clusters
        debug_report_global_energetics = self.debug_report_global_energetics
        debug_report_processes         = self.debug_report_processes
        debug_newton_method            = self.debug_newton_method
        debug_check_processes          = self.debug_check_processes
        debug_check_lattice            = self.debug_check_lattice

        f = open(filename, 'w+')

        col = 24
        print('', file = f)
        print('random_seed'.ljust(col), random_seed, file = f)
        if tol_dx_newton is not None:
            print('tol_dx_newton'.ljust(col), tol_dx_newton, file = f)
        if tol_rhs_newton is not None:
            print('tol_rhs_newton'.ljust(col), tol_rhs_newton, file = f)
        if max_newton_iter is not None:
            print('max_newton_iter'.ljust(col), max_newton_iter, file = f)
        if n_gauss_pts is not None:
            print('n_gauss_pts'.ljust(col), n_gauss_pts, file = f)
        print('', file = f)
        if ramp is not None:
            temperature = 'ramp  {0}  {1}'.format(temperature, ramp)
        print('temperature'.ljust(col), temperature, ' # [K]', file = f)
        print('pressure'.ljust(col), pressure, ' # [bar]', file = f)

        print('', file = f)
        print('n_gas_species'.ljust(col), len(gas_specs_names), file = f)
        gas_specs_string = 'gas_specs_names'.ljust(col)
        for item in gas_specs_names:
            gas_specs_string += (' {0:9s} '.format(item))
        print(gas_specs_string, file = f)

        energies_string = 'gas_energies'.ljust(col)
        for item in gas_energies:
            energies_string += (' {0:9.4f} '.format(item))
        print(energies_string + ' # [eV]', file = f)

        molec_weights_string = 'gas_molec_weights'.ljust(col)
        for item in gas_molec_weights:
            molec_weights_string += (' {0:9.4f} '.format(item))
        print(molec_weights_string + ' # [uma]', file = f)

        molar_fracs_string = 'gas_molar_fracs'.ljust(col)
        for item in gas_molar_fracs:
            molar_fracs_string += (' {0:9.4f} '.format(item))
        print(molar_fracs_string + ' # [-]', file = f)

        print('', file = f)
        print('n_surf_species'.ljust(col), len(surf_specs_names), file = f)

        surf_specs_string = 'surf_specs_names'.ljust(col)
        for item in surf_specs_names:
            surf_specs_string += (' {0:7s} '.format(item))
        print(surf_specs_string, file = f)

        if surf_specs_dent is None:
            surf_specs_dent = []
            for i in range(len(surf_specs_names)):
                surf_specs_dent.append(surf_specs_names[i].count('*'))

        specs_dent_string = 'surf_specs_dent'.ljust(col)
        for item in surf_specs_dent:
            specs_dent_string += (' {0:7d} '.format(item))
        print(specs_dent_string, file = f)

        print('', file = f)
        if isinstance(snapshots, tuple):
            snapshots = ' '.join(str(i) for i in snapshots)
        print('snapshots'.ljust(col), acquirement, snapshots, file = f)
        if isinstance(process_statistics, tuple):
            process_statistics = ' '.join(str(i) for i in process_statistics)
        print('process_statistics'.ljust(col), acquirement, process_statistics,
            file = f)
        if isinstance(species_numbers, tuple):
            species_numbers = ' '.join(str(i) for i in species_numbers)
        print('species_numbers'.ljust(col), acquirement, species_numbers,
            file = f)

        if energetics_lists:
            print('energetics_lists'.ljust(col), energetics_lists, file = f)
        if process_lists:
            print('process_lists'.ljust(col), process_lists, file = f)

        print('', file = f)
        if event_report is True:
            print('event_report'.ljust(col), 'on', file = f)
        elif event_report is False:
            print('event_report'.ljust(col), 'off', file = f)
        else:
            print('event_report'.ljust(col), event_report, file = f)

        print('', file = f)
        if max_steps is None:
            max_steps = 'infinity'
        print('max_steps'.ljust(col), max_steps, file = f)
        if max_time is None:
            max_time = 'infinity'
        print('max_time'.ljust(col), max_time, file = f)

        print('', file = f)
        print('wall_time'.ljust(col), wall_time, ' # [s]', file = f)

        print('', file = f)
        if restart is False:
            print('no_restart', file = f)

        if enable_stiffness_scaling:
            print('', file = f)
            print('enable_stiffness_scaling', file = f)
            if stiffness_scale_all:
                print('stiffness_scale_all', file = f)
            if check_every:
                print('check_every', check_every, file = f)
            if min_separation:
                print('min_separation', min_separation, file = f)
            if max_separation:
                print('max_separation', max_separation, file = f)
            if max_qequil_separation:
                print('max_qequil_separation', max_qequil_separation, file = f)
            if tol_part_equil_ratio:
                print('tol_part_equil_ratio', tol_part_equil_ratio, file = f)
            if stiffn_coeff_threshold:
                print('stiffn_coeff_threshold', stiffn_coeff_threshold,
                    file = f)
            if scaling_factor:
                print('scaling_factor', scaling_factor, file = f)
            if max_qequil_separation:
                print('cache_product_clusters', cache_product_clusters,
                    file = f)

        if debug_report_global_energetics:
            print('debug_report_global_energetics', file = f)
        if debug_report_processes:
            print('debug_report_processes', file = f)
        if debug_newton_method:
            print('debug_newton_method', file = f)
        if debug_check_processes:
            print('debug_check_processes', file = f)
        if debug_check_lattice:
            print('debug_check_lattice', file = f)

        print('', file = f)
        print('finish', file = f)

        f.close()

    # -------------------------------------------------------------------
    #  PRINT ENERGETICS INPUT
    # -------------------------------------------------------------------

    def print_energetics_input(self, filename = None):

        try: os.mkdir(self.dirname)
        except: pass

        if filename is None:
            filename = os.path.join(self.dirname, 'energetics_input.dat')

        clusters = self.clusters

        f = open(filename, 'w+')

        print('', file = f)
        print('energetics', file = f)
        print('', file = f)

        for i, cluster in enumerate(clusters):
            print(self.separation_string, file = f)
            print('', file = f)
            block = cluster.create_cluster_block()
            print(block, file = f)

        print(self.separation_string, file = f)
        print('', file = f)
        print('end_energetics', file = f)

        f.close()

    # -------------------------------------------------------------------
    #  PRINT MECHANISM INPUT
    # -------------------------------------------------------------------

    def print_mechanism_input(self, filename = None):

        try: os.mkdir(self.dirname)
        except: pass

        if filename is None:
            filename = os.path.join(self.dirname, 'mechanism_input.dat')

        steps = self.steps

        f = open(filename, 'w+')

        print('', file = f)
        print('mechanism', file = f)
        print('', file = f)

        for i, step in enumerate(steps):
            print(self.separation_string, file = f)
            print('', file = f)
            block = step.create_step_block()
            print(block, file = f)

        print(self.separation_string, file = f)
        print('', file = f)
        print('end_mechanism', file = f)
            
        f.close()

    # -------------------------------------------------------------------
    #  PRINT LATTICE INPUT
    # -------------------------------------------------------------------

    def print_lattice_input(self, filename = None):

        try: os.mkdir(self.dirname)
        except: pass

        if filename is None:
            filename = os.path.join(self.dirname, 'lattice_input.dat')

        active_sites = self.active_sites

        f = open(filename, 'w+')

        names = []
        max_coord = 0
        for a in active_sites:
            if a.name not in names:
                names += [a.name]
            if len(a.neighbors) > max_coord:
                max_coord = len(a.neighbors)

        names = sorted(names)

        col = 18
        print('', file = f)
        print('lattice explicit', file = f)
        print('', file = f)
        print('n_sites'.ljust(col), len(active_sites), file = f)
        print('max_coord'.ljust(col), max_coord, file = f)
        print('n_site_types'.ljust(col), len(names), file = f)

        names_string = ''
        for name in names:
            names_string += ('{}  '.format(name))
        print('site_type_names'.ljust(col), names_string, file = f)

        print('', file = f)
        print('lattice_structure', file = f)
        for a in active_sites:
            x, y = a.position_xy
            string = '{0:4d} {1:7.2f} {2:7.2f}'.format(a.index + 1, x, y)
            string += '   {0:15s} {1:3d}'.format(a.name, len(a.neighbors))
            for n in a.neighbors:
                string += '{0:6d}'.format(n+1)
            print(string, file = f)
        print('end_lattice_structure', file = f)

        print('', file = f)
        print('end_lattice', file = f)
            
        f.close()

################################################################################
# CLUSTER VARIANT
################################################################################

class ClusterVariant():

    def __init__(self,
                 name               = None,
                 site_types         = None,
                 graph_multiplicity = None,
                 cluster_eng        = None,
                 angles             = None,
                 no_mirror_images   = None,
                 absl_orientation   = None):

        self.name               = name
        self.site_types         = site_types
        self.cluster_eng        = cluster_eng
        self.graph_multiplicity = graph_multiplicity
        self.angles             = angles
        self.no_mirror_images   = no_mirror_images
        self.absl_orientation   = absl_orientation

################################################################################
# CLUSTER
################################################################################

class Cluster():

    def __init__(self,
                 name               = None,
                 configuration      = None,
                 neighboring        = None,
                 site_types         = None,
                 graph_multiplicity = None,
                 cluster_eng        = None,
                 angles             = None,
                 no_mirror_images   = False,
                 absl_orientation   = None,
                 variants           = None):

        if variants is None:
            variants = []

        self.name               = name
        self.configuration      = configuration
        self.neighboring        = neighboring
        self.site_types         = site_types
        self.graph_multiplicity = graph_multiplicity
        self.cluster_eng        = cluster_eng
        self.angles             = angles
        self.no_mirror_images   = no_mirror_images
        self.absl_orientation   = absl_orientation
        self.variants           = variants

    # -------------------------------------------------------------------
    #  ADD VARIANT
    # -------------------------------------------------------------------

    def add_variant(self, sites = None, eng = None, **kwargs):

        var = ClusterVariant(**kwargs)

        if sites is not None:
            var.site_types = sites
        if eng is not None:
            var.cluster_eng = eng

        self.variants.append(var)

    # -------------------------------------------------------------------
    #  CREATE CLUSTER BLOCK
    # -------------------------------------------------------------------

    def create_cluster_block(self):

        name               = self.name
        configuration      = self.configuration
        neighboring        = self.neighboring
        site_types         = self.site_types
        graph_multiplicity = self.graph_multiplicity
        cluster_eng        = self.cluster_eng
        angles             = self.angles
        no_mirror_images   = self.no_mirror_images
        absl_orientation   = self.absl_orientation
        variants           = self.variants
        block              = ''

        if type(site_types) is str:
            site_types = [site_types]

        if name is None:
            name = configuration.replace(' ', '')
        block += 'cluster  {}\n\n'.format(name)

        configuration = configuration.replace(' ', '')
        adsorbates = configuration.split('+')

        ads_config = ''
        site_types_config = []
        for i, ads in enumerate([ ads for ads in adsorbates if '*' in ads \
                                                            or '&' in ads ]):
            adsorbate = ads.split('*')[0]
            site = ads.split('*')[ -len(ads.split('*')) + 1 : ]
            if '&' in adsorbate:
                site = site[0]
                if variants == [] and site in ('&', ''):
                    raise TypeError('Specify the site type for &!')
                site_types_config.append(str(site))
                ads_config += '    & &        &\n'
            else:
                adsorbate += '*' * len(site)
                for j, site in enumerate(site):
                    site_types_config.append(site)
                    ads_config += '  {0:3d} {1:6s} {2:3d}\n'.format(i+1,
                        adsorbate, j+1)

        if site_types is None:
            site_types = site_types_config
        else:
            if len(site_types) != len(site_types_config):
                raise ValueError('Wrong number of site in site_types!')

        if site_types and len(site_types) != 1:
            neighboring = calculate_neighboring(neighboring, site_types)

        if site_types and not variants and graph_multiplicity is None:
            graph_multiplicity = calculate_multiplicity(adsorbates, site_types)

        block += '  sites {}\n'.format(len(site_types))
        if neighboring is not None:
            block += '  neighboring {}\n'.format(neighboring)
        block += '  lattice_state\n' + ads_config

        if variants:

            for var in variants:

                if type(var.site_types) is str:
                    var.site_types = [var.site_types]

                if var.name is None and var.site_types:
                    var.name = '_'.join(var.site_types)

                if len(site_types) != len(var.site_types):
                    raise ValueError('Wrong number of site in site_types!')

                var_site_types_string = ''
                for item in var.site_types:
                    var_site_types_string += (' {0:6s}'.format(item))

                if var.site_types and var.graph_multiplicity is None:
                    var.graph_multiplicity = calculate_multiplicity(adsorbates,
                                                                var.site_types)

                block += '\n  variant  {}\n'.format(var.name)
                block += '    site_types  {}\n'.format(var_site_types_string)
                if var.graph_multiplicity is not None:
                    block += '    graph_multiplicity  {}\n'.format(\
                        var.graph_multiplicity)
                block += '    cluster_eng {0:13.5e}\n'.format(var.cluster_eng)
                if var.angles is not None:
                    block += '    angles {}\n'.format(var.angles)
                if var.no_mirror_images is True:
                    block += '    no_mirror_images\n'
                if var.absl_orientation is not None:
                    block += '    absl_orientation {}\n'.format(\
                        var.absl_orientation)

                block += '  end_variant\n'

        else:
            site_types_string = ''
            for item in site_types:
                site_types_string += (' {0:6s}'.format(item))

            block += '  site_types  {}\n'.format(site_types_string)
            if graph_multiplicity is not None:
                block += '  graph_multiplicity  {}\n'.format(graph_multiplicity)
            block += '  cluster_eng {0:13.5e}\n'.format(cluster_eng)
            if angles is not None:
                block += '  angles {}\n'.format(angles)
            if no_mirror_images is True:
                block += '  no_mirror_images\n'
            if absl_orientation is not None:
                block += '  absl_orientation {}\n'.format(absl_orientation)

        block += '\nend_cluster\n'

        return block

################################################################################
# STEP VARIANT
################################################################################

class StepVariant():

    def __init__(self,
                 name              = None,
                 site_types        = None,
                 pre_expon         = None,
                 pre_expon_rev     = None,
                 pe_ratio          = None,
                 activ_eng         = None,
                 angles            = None,
                 no_mirror_images  = None,
                 absl_orientation  = None,):

        self.name             = name
        self.site_types       = site_types
        self.pre_expon        = pre_expon
        self.pre_expon_rev    = pre_expon_rev
        self.pe_ratio         = pe_ratio
        self.activ_eng        = activ_eng
        self.angles           = angles
        self.no_mirror_images = no_mirror_images
        self.absl_orientation = absl_orientation

################################################################################
# STEP
################################################################################

class Step():

    def __init__(self,
                 name               = None,
                 reaction           = None,
                 reversible         = False,
                 gas_reacs_prods    = None,
                 neighboring        = None,
                 site_types         = None,
                 pre_expon          = None,
                 pre_expon_rev      = None,
                 pe_ratio           = None,
                 activ_eng          = None,
                 prox_factor        = None,
                 angles             = None,
                 no_mirror_images   = None,
                 absl_orientation   = None,
                 stiffness_scalable = None,
                 variants           = None):

        if variants is None:
            variants = []

        self.name               = name
        self.reaction           = reaction
        self.reversible         = reversible
        self.gas_reacs_prods    = gas_reacs_prods
        self.neighboring        = neighboring
        self.site_types         = site_types
        self.pre_expon          = pre_expon
        self.pre_expon_rev      = pre_expon_rev
        self.pe_ratio           = pe_ratio
        self.activ_eng          = activ_eng
        self.prox_factor        = prox_factor
        self.angles             = angles
        self.no_mirror_images   = no_mirror_images
        self.absl_orientation   = absl_orientation
        self.stiffness_scalable = stiffness_scalable
        self.variants           = variants

    # -------------------------------------------------------------------
    #  ADD VARIANT
    # -------------------------------------------------------------------

    def add_variant(self, sites = None, A_fwd = None, A_rev = None,
                    Eact = None, **kwargs):

        var = StepVariant(**kwargs)

        if sites is not None:
            var.site_types = sites
        if A_fwd is not None:
            var.pre_expon = A_fwd
        if A_rev is not None:
            var.pre_expon_rev = A_rev
        if Eact is not None:
            var.activ_eng = Eact

        self.variants.append(var)

    # -------------------------------------------------------------------
    #  CREATE STEP BLOCK
    # -------------------------------------------------------------------

    def create_step_block(self):

        name               = self.name
        reaction           = self.reaction
        reversible         = self.reversible
        gas_reacs_prods    = self.gas_reacs_prods
        neighboring        = self.neighboring
        site_types         = self.site_types
        pre_expon          = self.pre_expon
        pre_expon_rev      = self.pre_expon_rev
        pe_ratio           = self.pe_ratio
        activ_eng          = self.activ_eng
        prox_factor        = self.prox_factor
        angles             = self.angles
        no_mirror_images   = self.no_mirror_images
        absl_orientation   = self.absl_orientation
        stiffness_scalable = self.stiffness_scalable
        variants           = self.variants
        block              = ''

        if name is None:
            name = reaction.replace(' ', '')

        reaction = reaction.replace(' ', '')

        if '<->' in reaction or '<=>' in reaction:
            reversible = True
            reactants_str, products_str = re.split('<->|<=>', reaction)
        elif '->' in reaction or '=>' in reaction:
            reactants_str, products_str = re.split('->|=>', reaction)
        else:
            raise NameError('no arrow in reaction string')

        if reversible is True:
            block += 'reversible_step {}\n\n'.format(name)
        else:
            block += 'step {}\n\n'.format(name)

        reactants = reactants_str.split('+')
        products = products_str.split('+')

        gas_reactants = []
        for i, react in enumerate([ r for r in reactants if '*' not in r ]):
            gas_reactants.append(react)

        gas_products = []
        for i, prod in enumerate([ p for p in products if '*' not in p ]):
            gas_products.append(prod)

        gas_react_string = '  gas_reacs_prods  '
        try:
            for item in gas_reactants:
                gas_react_string += ('{0:6s} {1:3d}'.format(item, -1))
        except:
            pass
        try:
            for item in gas_products:
                gas_react_string += ('{0:6s} {1:3d}'.format(item, 1))
        except:
            pass

        block += gas_react_string + '\n'

        react_conf = ''
        sites_react = []
        for i, react in enumerate([ r for r in reactants if '*' in r ]):
            sites = react.split('*')[1:]
            ads = react.split('*')[0] + '*' * len(sites)
            for j, site in enumerate(sites):
                sites_react.append(site)
                react_conf += '  {0:3d} {1:6s} {2:3d}\n'.format(i+1, ads, j+1)

        prod_conf = ''
        sites_prod = []
        for i, prod in enumerate([ p for p in products if '*' in p ]):
            sites = prod.split('*')[1:]
            ads = prod.split('*')[0] + '*' * len(sites)
            for j, site in enumerate(sites):
                sites_prod.append(site)
                prod_conf += '  {0:3d} {1:6s} {2:3d}\n'.format(i+1, ads, j+1)

        if sites_react != sites_prod:
            raise ValueError('Different sites in initial and final configs!')
        if site_types is None:
            site_types = sites_react
        else:
            if len(site_types) != len(sites_react):
                raise ValueError('Different number of sites in step!')

        if site_types and len(site_types) != 1:
            neighboring = calculate_neighboring(neighboring, site_types)

        block += '  sites {}\n'.format(len(site_types))
        if neighboring is not None:
            block += '  neighboring {}\n'.format(neighboring)
            
        block += '  initial\n' + react_conf + '  final\n' + prod_conf

        if variants:

            for var in variants:

                if type(var.site_types) is str:
                    var.site_types = [var.site_types]

                if var.name is None and var.site_types:
                    var.name = '_'.join(var.site_types)

                if len(site_types) != len(var.site_types):
                    raise ValueError('Wrong number of site in site_types!')

                var_site_types_string = ''
                for item in var.site_types:
                    var_site_types_string += (' {0:6s}'.format(item))

                block += '\n  variant  {}\n'.format(var.name)
                block += '    site_types  {}\n'.format(var_site_types_string)

                if var.pre_expon_rev is not None:
                    var.pe_ratio = var.pre_expon / var.pre_expon_rev

                if type(var.pre_expon) in (int, float):
                    block += '    pre_expon  {0:13.5e}'.format(var.pre_expon)
                elif len(var.pre_expon) == 7:
                    block += '    pre_expon  '
                    for n in var.pre_expon:
                        if n is None:
                            block += ' 0.'
                        else:
                            block += '{0:13.5e}'.format(n)
                else:
                    raise ValueError( 'Wrong length of pre exponential!')
                block += '\n'
                if reversible is True:
                    block += '    pe_ratio   {0:13.5e}\n'.format(var.pe_ratio)
                block += '    activ_eng  {0:13.5e}\n'.format(var.activ_eng)

                if var.angles is not None:
                    block += '  angles {}\n'.format(var.angles)
                if var.no_mirror_images is not None:
                    block += '  no_mirror_images\n'
                if var.absl_orientation is not None:
                    block += '  absl_orientation {}\n'.format( \
                        var.absl_orientation)
                
                block += '  end_variant\n'

        else:

            if pre_expon_rev is not None:
                pe_ratio = pre_expon / pre_expon_rev

            site_types_string = ''
            for item in site_types:
                site_types_string += (' {0:6s}'.format(item))
            block += '  site_types  {}\n'.format(site_types_string)
            if type(pre_expon) in (int, float):
                block += '  pre_expon  {0:13.5e}'.format(pre_expon)
            elif len(pre_expon) == 7:
                block += '  pre_expon  '
                for n in pre_expon:
                    if n is None:
                        block += ' 0.'
                    else:
                        block += '{0:13.5e}'.format(n)
            else:
                raise ValueError('Wrong length of pre exponential!')
            block += '\n'
            if reversible is True:
                block += '  pe_ratio   {0:13.5e}\n'.format(pe_ratio)
            block += '  activ_eng  {0:13.5e}\n'.format(activ_eng)

            if angles is not None:
                block += '  angles {}\n'.format(angles)
            if no_mirror_images is not None:
                block += '  no_mirror_images\n'
            if absl_orientation is not None:
                block += '  absl_orientation {}\n'.format(absl_orientation)

        if prox_factor is not None:
            block += '  prox_factor {}\n'.format(prox_factor)
        if stiffness_scalable is not None:
            block += '  stiffness_scalable {}\n'.format(stiffness_scalable)

        if reversible is True:
            block += '\nend_reversible_step\n'
        else:
            block += '\nend_step\n'

        return block

################################################################################
# CALCULATE NEIGHBORING
################################################################################

def calculate_neighboring(neighboring, site_types):

    if neighboring in (None, 'all'):
        vector_str = [ str(i) for i in range(1, len(site_types)+1) ]
        neigh = list(itertools.combinations(vector_str, 2))
        for i in range(len(neigh)):
            neigh[i] = '-'.join(neigh[i])
        neighboring = ' '.join(neigh)
    elif neighboring == 'line':
        neigh = [ str(i) + '-' + str(i+1) for i in \
           range(1, len(site_types)) ]
        neighboring = ' '.join(neigh)
    elif neighboring == 'circle':
        neigh = [ str(i) + '-' + str(i+1) for i in \
           range(1, len(site_types)) ]
        neighboring = ' '.join(neigh) + ' {}-1'.format(len(site_types))

    return neighboring

################################################################################
# CALCULATE MULTIPLICITY
################################################################################

def calculate_multiplicity(adsorbates, site_types):

    graph_multiplicity = None

    if len(adsorbates) == 2 and adsorbates[0] == adsorbates[1]:
        half_len = int(len(site_types)/2)
        if list(site_types[:half_len]) == list(site_types[half_len:]) \
        or list(site_types[:half_len]) == list(site_types[half_len:])[::-1]:
            graph_multiplicity = 2

    #TODO: EXTEND TO OTHER CASES

    return graph_multiplicity

################################################################################
# END
################################################################################
