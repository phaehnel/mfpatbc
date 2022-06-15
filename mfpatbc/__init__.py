#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The mfpatbc module allows to create input data for a phase-averaged tidal 
boundary condition as required by the Python package FloPy.
"""

from .phase_average_tbc import *
from .modflow_bc import *

# Load regression model coefficients
import pandas as pd

path_module = __file__[:-11]


major = 1
minor = 0
micro = 0
__version__ = f'{major}.{minor}.{micro}'

__pakname__ = 'mfpatbc'

__author__ = {'Patrick Haehnel': 'patrick.haehnel@uni-oldenburg.de'}