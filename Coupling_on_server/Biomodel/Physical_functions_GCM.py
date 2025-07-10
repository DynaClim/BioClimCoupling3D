#!/usr/bin/env python
# coding: utf-8

import numpy as np
from Constants_GCM import *

def atm_matter_quantity(P_profile, T_profile,Z,surface, g=3.71):
    """
    Calculates the total amount of matter (in moles) of a chemical species in the atmosphere.

    Parameters:
    - P_profile : array, partial pressure (Pa) of the species (from top to surface)
    - T_profile : array, temperature (K), in the same order as P_profile
        - Z : array, altitudes (m)
    - surface : surface area (in m²)
    - g : gravity (m/s²), default is 3.71 m/s² (Mars)

    Returns:
    - n_total : total amount of matter (mol) of the species
    """
    dZ = np.abs(np.diff(Z))
    # Calculate the mean temperatures and pressures between layers
    T_mean = (T_profile[:-1] + T_profile[1:]) / 2
    P_mean = (P_profile[:-1] + P_profile[1:]) / 2
    
    # Calculate the volume of each atmospheric layer (cylindrical volume)
    V_shell = dZ * surface

    # Calculate the amount of matter in each layer
    n_shell = (P_mean * V_shell) / (R * T_mean)

    # Sum the total amount of matter across all layers
    n_total = np.sum(n_shell)

    return n_total

