#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Patrick Haehnel
2022-01-29

Contains functions to calculate hydraulic heads for phase-averaged tidal boundary 
condition at the aquifer-ocean interface in the intertidal zone.

Based on an analytical solution presented by Nuttle (1991), hlim, which describes
phase-averaged intertidal head distribution under tidal forcing for very shallow 
sloping beaches (pavg_heads_lim).

Empirical correction function adapts this solutions to deviations from the
assumptions of this solution due to parameters 
    - horizontal hydraulic conductivity (hk) [m/s]
    - specific yield (sy) [-]
    - saturated aquifer thickness (D) [m]
        D = hs - z_D, hs -> phase-averaged sea level [m asl], 
        z_D -> aquifer base [m asl]
    - vertical anisotropy (vka = vk/hk, vk is vertical hydraulic conductivity) 
    - beach slope (beta)
using a linear regression model describing the phase-averaged exit point 
elevation of the groundwater table (zep). This is the elevation at which
groundwater table is equal to surface elevation on phase-averaged basis
(predict_zep).

Adaption of hlim is a cutoff threshold where zep is used instead of hlim
whenever hlim > zep (pavg_heads).

Additionally, a function is included allowing a user defined, constant setting of
the empirical correction function (pavg_heads_zep_const).


References:
    Nuttle, W. K. (1991). Comment on “Tidal dynamics of the water table in 
        beaches” by Peter Nielsen. Water Resources Research, 27(7), 1781–1782. 
        https://doi.org/10.1029/91WR00939
"""

import numpy as np


def pavg_heads_lim(hs, A, z):
    """
    Phase-averaged intertidal hydraulic heads according to the model of
    Nuttle (1991).
    
    Additionally sets heads beloe low water to hs and heads above high water
    to high water (hs+A).

    Parameters
    ----------
    hs : float
        Phase-averaged sea level.
    A : float
        Tidal amplitude.
    z : float or array of floats
        Surface elevation.

    Returns
    -------
    havg : float or array of floats
        Phase-averaged hydraulic heads.

    """
    
    # Normalized elevation, relative to MSL and amplitude
    z_norm = (z-hs)/A
    acosz = np.arccos(z_norm)
    
    # Calculate the limiting case intertidal phase-average heads
    havg = 1/np.pi * ((hs*acosz) + z*(np.pi - acosz) +\
                      A*np.sqrt(1 - z_norm**2))
        
    # Set heads below low water to hs (use <= and >=, otherwise nan can occur
    # in the intertidal zone if z very close to hs-A or hs+A)
    havg[z <= (hs - A)] = hs
    havg[z >= (hs + A)] = hs + A
        
    return havg


def pavg_heads_zep_const(hs, A, z, zep_norm_const):
    """
    Calculate phase-averaged intertidal hydraulic heads with a constant user
    defined value for the empirical correction function.

    Parameters
    ----------
    hs : float
        Phase-averaged sea level.
    A : float
        Tidal amplitude.
    z : float or array of floats
        Surface elevation.
    zep_norm_const : float
        User defined value for empirical correction function. Must be in [0, 1].

    Returns
    -------
    havg : float or array of floats
        Phase-averaged hydraulic heads.

    """
    
    # Get limiting case heads for the intertidal
    havg = pavg_heads_lim(hs, A, z)
    
    # Transform relative to amplitude
    zep = hs + (A * zep_norm_const)
    
    # Set heads above exit point to exit point
    above_zep = havg > zep
    havg[above_zep] = zep
    
    return havg
    

def predict_zep(hs, A, T, hk, sy, D, vka, slope, lm):
    """
    Calculate value(s) of empirical correction function using the linear 
    regression model and derives phase-averaged exit point elevation from this.

    Parameters
    ----------
    hs : float
        Phase-averaged sea level.
    A : float
        Tidal amplitude.
    T : float
        Period length of tidal cycle.
    hk : float
        Horizontal hydraulic conductivity, dimensional.
    sy : float
        Specific yield.
    D : float
        Depth of the aquifer base below mean sea level (z_D in paper!), dimensional.
    vka : float
        vertical anisotropy.
    slope : float
        beach slope as gradient, tan(beta).
    lm : pandas DataFrame
        Terms and coefficients of the linear regression model.
        Terms are in R formula style and changed to work in Python

    Returns
    -------
    zep : float
        Estimated phase-averaged exit point elevation, cutoff threshhold.

    """
    
    # Convert parameters for empirical model
    hk_norm = hk * T / A
    D_norm = D / A
    
    # Extract variables and adapt for use in Python, all except intercept
    # Turn interaction signs to multiplication
    # Adapt power sign
    # delete I from e.g. polynomial terms
    terms = lm.iloc[::, 0].values[1::]
    terms = [tk.replace(':', '*') for tk in terms]
    terms = [tk.replace('^', '**') for tk in terms]
    terms = [tk.replace('I', '') for tk in terms]
    
    # All coeffiecients but intercept
    coef = lm.Coefficients.values[1::]
    
    log = np.log
    tan = np.tan
    
    beta_rad = np.arctan(slope)
    
    # Calculate phase-average exit point from empirical model, normalized 
    # Extract intercept seperatly
    intercept = lm.Coefficients.values[0]
    
    # Calculate remaining terms using the input file information
    scope = locals()
    zep_norm_transform = intercept
    for (c, tk) in zip(coef, terms):
        zep_norm_transform += c * eval(tk, scope)
    
    # Retransform using logistic function
    zep_norm = 1 / (1 + np.exp(-zep_norm_transform))
    
    # Transform relative to amplitude
    zep = hs + (A * zep_norm)
    
    return zep
    
    
def pavg_heads(hs, A, z, T, hk, sy, D, vka, slope, lm):
    """
    Calculates phase-averaged intertidal hydraulic heads for phase-averaged
    tidal boundary condition (PA-TBC) using empirical correction function.

    Parameters
    ----------
    hs : float
        Phase-averaged sea level.
    A : float
        Tidal amplitude.
    z : float or array of floats
        Surface elevation.
    T : float
        Period length of tidal cycle.
    hk : float
        Horizontal hydraulic conductivity, dimensional.
    sy : float
        Specific yield.
    D : float
        Depth of the aquifer base below mean sea level (z_D in paper!), dimensional.
    vka : float
        vertical anisotropy.
    slope : float
        beach slope as gradient, tan(beta).
    lm : pandas DataFrame
        Terms and coefficients of the linear regression model.

    Returns
    -------
    havg : float or array of floats
        Phase-averaged hydraulic heads.

    """
    
    # Get limiting case heads for the intertidal
    havg = pavg_heads_lim(hs, A, z)
    
    # Get prediction for phase-averaged exit point elevation
    zep = predict_zep(hs, A, T, hk, sy, D, vka, slope, lm)
    
    # Set heads above exit point to exit point
    above_zep = havg > zep
    # Distinguish between cases with all single parameter values and
    # such with at least one spatially changing parameter
    havg[above_zep] = zep if zep.size == 1 else zep[above_zep]
    
    return havg

    
    