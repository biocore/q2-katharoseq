# ----------------------------------------------------------------------------
# Copyright (c) 2022--, katharoseq development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from . import _version
from ._methods import read_count_threshold, estimating_biomass, control_type
from ._methods import biomass_plot

__version__ = _version.get_versions()['version']
__all__ = ['read_count_threshold', 'estimating_biomass', 'biomass_plot',
           'control_type']
