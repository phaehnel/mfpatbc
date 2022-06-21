#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Patrick Haehnel
2022-05-17

Contains functions to calculate hydraulic heads for phase-averaged tidal boundary 
condition at the aquifer-ocean interface in the intertidal zone.

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
    hs : float or array of floats (nrow, ncol)
        Phase-averaged sea level.
    A : float or array of floats (nrow, ncol)
        Tidal amplitude.
    z : float or array of floats (nrow, ncol)
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
    below_lw = z <= (hs - A)
    above_hw = z >= (hs + A)
    havg[below_lw] = hs if hs.size == 1 else hs[below_lw]
    havg[above_hw] = hs + A if hs.size == 1 else hs[above_hw] + A[above_hw]
        
    return havg


def pavg_heads_zep_const(hs, A, z, zep_norm_const):
    """
    Calculate phase-averaged intertidal hydraulic heads with a constant user
    defined value for the empirical correction function.

    Parameters
    ----------
    hs : float or array of floats (nrow, ncol)
        Phase-averaged sea level.
    A : float or array of floats (nrow, ncol)
        Tidal amplitude.
    z : float or array of floats (nrow, ncol)
        Surface elevation.
    zep_norm_const : float or array of floats (nrow, ncol)
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
    havg[above_zep] = zep if zep.size == 1 else zep[above_zep]
    
    return havg
    

def predict_zep(hs, A, T, hk, sy, D, vka, slope, lm):
    """
    Calculate value(s) of empirical correction function using the linear 
    regression model and derives phase-averaged exit point elevation from this.

    Parameters
    ----------
    hs : float or array of floats (nrow, ncol)
        Phase-averaged sea level.
    A : float or array of floats (nrow, ncol)
        Tidal amplitude.
    T : float
        Period length of tidal cycle.
    hk : float or array of floats (nrow, ncol)
        Horizontal hydraulic conductivity, dimensional [m/s].
    sy : float or array of floats (nrow, ncol)
        Specific yield.
    D : float or array of floats (nrow, ncol)
        Saturated aquifer thickness (distance phase-averaged sea level to 
        elevation of aquifer base), dimensional [m].    
    vka : float or array of floats (nrow, ncol)
        vertical anisotropy.
    slope : float or array of floats (nrow, ncol)
        beach slope as gradient, tan(beta).
    lm : pandas DataFrame
        Terms and coefficients of the linear regression model.
        Terms are in R formula style and changed to work in Python

    Returns
    -------
    zep : float or array of floats (nrow, ncol)
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
    hs : float or array of floats (nrow, ncol)
        Phase-averaged sea level.
    A : float or array of floats (nrow, ncol)
        Tidal amplitude.
    z : float or array of floats (nrow, ncol)
        Surface elevation.
    T : float or array of floats (nrow, ncol)
        Period length of tidal cycle.
    hk : float or array of floats (nrow, ncol)
        Horizontal hydraulic conductivity, dimensional [m/s].
    sy : float or array of floats (nrow, ncol)
        Specific yield.
    D : float or array of floats (nrow, ncol)
        Saturated aquifer thickness (distance phase-averaged sea level to 
        elevation of aquifer base), dimensional [m].
    vka : float or array of floats (nrow, ncol)
        vertical anisotropy.
    slope : float or array of floats (nrow, ncol)
        beach slope as gradient, tan(beta).
    lm : pandas DataFrame
        Terms and coefficients of the linear regression model.

    Returns
    -------
    havg : float or array of floats (nrow, ncol)
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

    
    