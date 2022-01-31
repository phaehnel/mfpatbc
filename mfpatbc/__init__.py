#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Patrick Hähnel
2022-01-29
"""

from .phase_average_tbc import *
from .modflow_bc import *

# Load regression model coefficients
import pandas as pd

path_module = __file__[:-11]