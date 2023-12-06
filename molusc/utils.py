from datetime import datetime as dt
import datetime
import warnings
import os, pathlib, sys
import logging

import numpy as np
import scipy as scipy
import scipy.stats as stats
import astropy.units as u
from astropy.table import Table
from time import time
from astropy.utils.exceptions import AstropyWarning
from multiprocessing import Process
from multiprocessing.pool import Pool
import multiprocessing as mp

# from pkgcore.config import load_config
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

def get_pro_sep(T_init, period, phase, ecc, a_peri, cos_i, semi_maj_a):
    """
    Calculate the projected separation given input orbital parameters

    """


    # Calculate projected separation for each generated companion:
    # 1. Calculate mean anomaly
    M = 2 * np.pi * T_init / period - phase
    # 2. Calculate eccentric anomaly iteratively
    prev_E = 0.0
    current_E = M
    while abs(current_E - prev_E) > 0.00001:
        prev_E = current_E
        current_E = M + ecc * np.sin(prev_E)
    
    # 3. Calculate true anomaly
    f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - ecc) / (1 + ecc)))
    
    # 4. Calculate projected separation in AU
    alpha = f + a_peri
    sqt = np.sqrt(np.sin(alpha)**2+np.cos(alpha)**2 * cos_i**2)
    pro_sep = semi_maj_a * (1-ecc**2)/(1+ecc*np.cos(f))*sqt
    
    return pro_sep

def calc_anomaly_gaia(num_generated,T_0,period,phase,ecc,arg_peri,cos_i,semi_maj_a):
    # Calculate projected separation for each generated companion
    f = np.zeros(num_generated)
    f_ecc_term = np.sqrt((1 - ecc) / (1 + ecc))
    pi2per = 2*np.pi/period

    for i in range(num_generated):
        # TODO: This should care about the date we observed the system
        # 1. Calculate mean anomaly
        M = pi2per[i] * T_0 - phase[i]
        # 2. Calculate eccentric anomaly iteratively
        prev_E = 0.0
        current_E = M
        while abs(current_E - prev_E) > 0.00001:
            prev_E = current_E
            current_E = M + ecc[i] * np.sin(prev_E)
        # 3. Calculate true anomaly
        f[i] = 2 * np.arctan2(np.tan(current_E / 2), f_ecc_term[i])
    # 4. Calculate projected separation in AU
    alpha = f + arg_peri
    sqt = np.sqrt(np.sin(alpha)**2+np.cos(alpha)**2 * cos_i**2)
    pro_sep = semi_maj_a * (1-ecc**2)/(1+ecc*np.cos(f))*sqt

    return pro_sep

def calculate_RV(period, mass_ratio, a, ecc, cos_i, arg_peri, phase, MJD):
    # Calculates the RVs for each item when passed arrays of orbital parameters
    # Inputs: Arrays of Period, Mass Ratio, Semi-Major Axis, eccentricity, inclination, arg peri, phase, calculation times
    # Outputs: Velocity Semi-Amplitude (km/s), RVs at each time in MJD
    
    sin_i = np.sin(np.arccos(cos_i))

    n = len(period)
    ndates = len(MJD)
    # Create a blank array with one row for every input period
    # and each row contains RVs equivalent to the observation dates
    RV = np.zeros((n,ndates))
    
    a_star = a * (mass_ratio / (mass_ratio + 1))

    K = (2*np.pi/period) * (a_star * sin_i) / np.sqrt(1-ecc**2) # AU/days
    K = K * 1731.48  # km/s

    for i in range(n):  # Iterate over companions
        for j in range(0, len(MJD)):  # Iterate over times
            # Find E
            M = 2 * np.pi * MJD[j] / period[i] - phase[i]
            prev_E = 0.0
            current_E = M

            while abs(current_E - prev_E) > 0.00001:
                prev_E = current_E
                current_E = M + ecc[i] * np.sin(prev_E)

            # Find true anomaly, f
            f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - ecc[i]) / (1 + ecc[i])))
            # Find predicted RV
            RV[i][j] = K[i] * (np.sin(arg_peri[i] + f) + ecc[i] * np.sin(arg_peri[i]))  # km/s

    return K, RV

def calculate_unscaled_RV(period, ecc, arg_peri, phase, MJD):
    """
    Calculate the *unscaled* radial velocity for a given period,
    eccentricity, arg_periastron, and phase at the given MJD dates.
    Arrays are NOT accepted for any argument except MJD. 

    Returns only the unscaled RV. 
    """
    
    ndates = len(MJD)
    
    f_ecc_term = np.sqrt((1 - ecc) / (1 + ecc))
    pi2per = 2*np.pi/period

    true_anomaly = np.zeros(ndates)

    for j in range(0, len(MJD)):  # Iterate over times
        # Find E
        M = pi2per * MJD[j] - phase
        prev_E = 0.0
        current_E = M

        while abs(current_E - prev_E) > 0.00001:
            prev_E = current_E
            current_E = M + ecc * np.sin(prev_E)

        # Find true anomaly, f
        true_anomaly[j] = 2 * np.arctan2(np.tan(current_E / 2), f_ecc_term)
        # Find predicted RV
    RV = np.sin(arg_peri + true_anomaly) + ecc * np.sin(arg_peri)  # km/s

    return RV
