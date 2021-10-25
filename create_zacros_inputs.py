#!/usr/bin/env python

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import numpy as np
from zacros_inputs import (ZacrosInputs, Cluster, Step, ClusterVariant,
                           StepVariant)

################################################################################
# ZACROS INPUTS
################################################################################

zin = ZacrosInputs(dirname = '.')

#-------------------------------------------------------------------------------
# SIMULATION INPUT
#-------------------------------------------------------------------------------

zin.random_seed        = 154686

zin.temperature        = 600.0 # [K]
zin.ramp               = None # [K/s]
zin.pressure           = 1.0 # [bar]

zin.gas_specs_names    = [    'A',     'B',     'C']
zin.gas_energies       = [  0.000,   0.000,   0.000] # [eV]
zin.gas_molec_weights  = [  1.000,   2.000,   3.000] # [uma]
zin.gas_molar_fracs    = [  0.500,   0.500,   0.000] # [-]

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

# Print simulation_input.dat
zin.print_simulation_input()

#-------------------------------------------------------------------------------
# ENERGETICS INPUT
#-------------------------------------------------------------------------------

# Binding Energies of A*
cls = Cluster(configuration = 'A*')
cls.add_variant(sites = 'top', eng = -1.0)
cls.add_variant(sites = 'brg', eng = -1.2)
zin.clusters.append(cls)

# Binding Energies of B*
cls = Cluster(configuration = 'B*')
cls.add_variant(sites = 'top', eng = -1.5)
cls.add_variant(sites = 'brg', eng = -1.7)
zin.clusters.append(cls)

# Lateral Interactions of A* and A*
cls = Cluster(configuration = 'A* + A*')
cls.add_variant(sites = ('top', 'top'), eng = +0.1)
cls.add_variant(sites = ('top', 'brg'), eng = +0.1)
cls.add_variant(sites = ('brg', 'brg'), eng = +0.1)
zin.clusters.append(cls)

# Lateral Interactions of B* and B*
cls = Cluster(configuration = 'B* + B*')
cls.add_variant(sites = ('top', 'top'), eng = +0.1)
cls.add_variant(sites = ('top', 'brg'), eng = +0.1)
cls.add_variant(sites = ('brg', 'brg'), eng = +0.1)
zin.clusters.append(cls)

# Lateral Interactions of A* and B*
cls = Cluster(configuration = 'A* + B*')
cls.add_variant(sites = ('top', 'top'), eng = +0.1)
cls.add_variant(sites = ('top', 'brg'), eng = +0.1)
cls.add_variant(sites = ('brg', 'brg'), eng = +0.1)
zin.clusters.append(cls)

# Print energetics_input.dat
zin.print_energetics_input()

#-------------------------------------------------------------------------------
# MECHANISM INPUT
#-------------------------------------------------------------------------------


# Adsorption of A
stp = Step(reaction = 'A + * <=> A*')
stp.add_variant(sites = 'top', k0 = 1., k0_rev = 1., Eact = 2.)
stp.add_variant(sites = 'brg', A_fwd = 1., A_rev = 1., Eact = 2.)
zin.steps.append(stp)

# Adsorption of A
stp = Step(reaction = 'A + * <=> A*')
stp.add_variant(sites = 'top', A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = 'brg', A_fwd = 1., A_rev = 1., Eact = 2.)
zin.steps.append(stp)

# Adsorption of B
stp = Step(reaction = 'B + * <=> B*')
stp.add_variant(sites = 'top', A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = 'brg', A_fwd = 1., A_rev = 1., Eact = 2.)
zin.steps.append(stp)

# Surface Diffusion of A*
stp = Step(reaction = 'A* + * <=> * + A*')
stp.add_variant(sites = ('top', 'top'), A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = ('top', 'brg'), A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = ('brg', 'top'), A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = ('brg', 'brg'), A_fwd = 1., A_rev = 1., Eact = 2.)
zin.steps.append(stp)

# Surface Diffusion of B*
stp = Step(reaction = 'B* + * <=> * + B*')
stp.add_variant(sites = ('top', 'top'), A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = ('top', 'brg'), A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = ('brg', 'top'), A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = ('brg', 'brg'), A_fwd = 1., A_rev = 1., Eact = 2.)
zin.steps.append(stp)

# Reaction between A* and B*
stp = Step(reaction = 'A* + B* <=> * + * + C')
stp.add_variant(sites = ('top', 'top'), A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = ('top', 'brg'), A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = ('brg', 'top'), A_fwd = 1., A_rev = 1., Eact = 2.)
stp.add_variant(sites = ('brg', 'brg'), A_fwd = 1., A_rev = 1., Eact = 2.)
zin.steps.append(stp)

# Print mechanism_input.dat
zin.print_mechanism_input()

################################################################################
# END
################################################################################
