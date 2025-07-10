#!/usr/bin/env python
# coding: utf-8

from Constants_GCM import *
import numpy as np


# ===========================================================
# Defining constants --> H2 , N2 and CH4
# ===========================================================

# Based on Kharecha et al, 2005

# Solubilities
alphaH     = 7.8e-4  # mol.L-1.bar-1
alphaG     = 1.4e-3  # mol.L-1.bar-1
alphaN     = 7.0e-4  # mol.L-1.bar-1
# Thermal diffusivities
thdH       = 5.0e-5  # cm2.s-1
thdG       = 1.8e-5  # cm2.s-1
#Piston velocities # Assuming no dependance on temperature
pvH        = 1.3e-2  # cm.s-1
pvG        = 4.5e-3  # cm.s-1
pvN        = 4.8e-3  # cm.s-1

# ===========================================================
# Defining constants --> CO2
# ===========================================================

def alphaC(T):
    """
    Computes the CO2 solubility depending on temperature

    Params:
     - T : float, temperature (K)

    Returns:
     - float, CO2 solubility (mol.L-1.bar-1)"""
    return np.exp(9345.17 / T - 167.8108 + 23.3585 * np.log(T) + (0.023517 - 2.3656e-4 * T + 4.7036e-7 * T**2) * 35.0)

pvC = 4.8e-3   # cm.s-1


def interface_flow(vp,alpha,p,c):
    """
    Computes the flux of a chemical specy at the interface between the ocean and the atmosphere
    Based on Kharecha et al, 2005

    Params:
     - vp : float, piston velocity (cm/s)
     - alpha : float, solubility (mol.L-1.bar-1)
     - p : float, partial pressure of the specy at the surface (bar)
     - c : float, concentration in the ocean (mol/L)

    Returns:
     - float, flow at the interface atmosphere -> ocean (mol/cm²/s)"""
    return vp*(alpha*p-c)*1e-3







