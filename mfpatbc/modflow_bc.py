#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Patrick Haehnel
2022-01-29
"""

import numpy as np
import pandas as pd
import mfpatbc


class PATBC():
    
    def __init__(self, m, hs, A, T, hk, vka, sy, z_D, slope, allow = None,
                 time_factor = 1, length_factor = 1):
        
        # Assign attributes
        self.nrow = m.nrow
        self.ncol = m.ncol
        self.z = m.dis.top.array
        self.hs = hs
        self.A = A
        self.T = T
        self.hk = hk
        self.vka = vka
        self.sy = sy
        self.z_D = z_D
        self.slope = slope
        self.allow = allow
        self.time_factor = time_factor
        self.length_factor = length_factor
        self.active = m.bas6.ibound.array[0]
        self.bottom = m.dis.botm.array[0]
        self.delr = m.dis.delr.array
        self.delc = m.dis.delc.array
        
        # If not provided all top layer cells will considered for PA-TBC
        if self.allow is None:
            self.allow = np.full((self.nrow, self.ncol), True)
            
        # Make properties spatial arrays if provided as scalar
        self.hk = self.make_array(self.hk)
        self.vka = self.make_array(self.vka)
        self.sy = self.make_array(self.sy)
        self.z_D = self.make_array(self.z_D)
        self.slope = self.make_array(self.slope)
        
        # Calculations required
        self.vertical_conductance = self.get_vertical_conductance()
    
    
    def get_drn_stress_period_data(self):
        
        # Find DRN cells
        is_drn = (self.hghb <= self.z) & self.allow & self.active
        
        # Get row, col and lay index of DRN cells
        row, col = np.where(is_drn)
        lay = np.full(row.shape, 0)
        
        # Select respective heads and conductances
        stage = self.z[row, col]
        cond = self.vertical_conductance[row, col]
        
        drn_stress_period_data = np.array((lay, row, col, stage, cond)).T
        
        return drn_stress_period_data
    
    
    def get_ghb_stress_period_data(self):
        
        lm = pd.read_csv(mfpatbc.path_module + 'coef_lm_sim-slopesteady_nwt.csv')
        
        # Saturated aquifer thickness
        D = self.hs - self.z_D
        
        # Calculate phase-averaged intertidal heads
        self.hghb = mfpatbc.pavg_heads(
            self.hs, self.A, self.z, self.T, self.hk, self.sy, D, self.vka, 
            self.slope, lm
            )
        
        # Find GHB cells
        is_ghb = (self.hghb > self.z) & self.allow & self.active
        
        # Get row, col and lay index of GHB cells
        row, col = np.where(is_ghb)
        lay = np.full(row.shape, 0)
        
        # Select respective heads and conductances
        stage = self.hghb[row, col]
        cond = self.vertical_conductance[row, col]
        
        ghb_stress_period_data = np.array((lay, row, col, stage, cond)).T
        
        return ghb_stress_period_data
            
        
    def get_vertical_conductance(self):
        
        length_flowpath = 0.5 * (self.z - self.bottom)
        
        DR, DC = np.meshgrid(self.delr, self.delc)
        area = DR * DC
        
        cond = (self.hk * area) / length_flowpath *\
            (self.time_factor / self.length_factor)
        
        return cond
    
    
    def make_array(self, parameter):
        if not hasattr(parameter, '__len__'):
            parameter = np.full((self.nrow, self.ncol), parameter)
            
        return parameter