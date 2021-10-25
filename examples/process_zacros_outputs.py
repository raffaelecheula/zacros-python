#!/usr/bin/env python

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
from zacros_outputs import ZacrosOutputs
import numpy as np

################################################################################
# ZACROS OUTPUTS
################################################################################

zout = ZacrosOutputs(dirname = '.')

zout.read_general_output()
zout.read_specnum_output()
zout.read_lattice_output()
zout.read_procstat_output()
zout.read_history_output()

filters = None

zout.plot_gas_prod_rates(filters = filters)
zout.plot_surf_coverages(filters = filters)
zout.plot_reaction_rates(filters = filters)

################################################################################
# END
################################################################################
