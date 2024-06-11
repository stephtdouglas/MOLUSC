import pathlib, os
import warnings

import h5py
import pytest
import numpy as np
import astropy.io.ascii as at
from astropy.table import Table
import scipy.stats as stats
import scipy

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
from molusc.ruwe import RUWE
from molusc.companions import Companions
from molusc.utils import set_null_limits

warnings.simplefilter('error', category=RuntimeWarning)

limits = set_null_limits()
pd_mu = 5.03
pd_sig = 2.28
q_exp = 0.0

# Use the demo star from the github page (HD 120411/HIP 67522)
# This is also the star shown in Figure 
star_ra = '13h50m06.28s '
star_dec = '-40d50m08.9s'
star_age = 0.015
star_mass = 1.3
plx = 8.017 #mas
star_ruwe = 0.91 #Rizzuto+ 2020

dist_pc = 1000/plx 
dist_au = dist_pc * 206265

# Compute the limits of the RUWE grid
log_sep_lims = np.array([1.275,3.325])
sep_lims = 10**log_sep_lims # mas
sep_lims_rad = np.radians(sep_lims/3.6e6)
# Small angle approximation means that we can multiply dist * radians
sep_lims_au = sep_lims_rad * dist_au

ref_file = os.path.join(repo_path,"reference_data/RuweTableGP.txt")
ref = at.read(ref_file)
nref = len(ref)

### This is Mackenna's original code, which uses interp2d ###
# Read in the reference table
file_name = ref_file
t = Table.read(file_name, format='ascii', delimiter=' ')
ruwe_dist = t

# f. Get predicted RUWE
#    convert from mas to AU
star_distance = 1 / (plx) * 2.063e+8 # AU
ruwe_dist['Sep(AU)'] = [star_distance * np.tan(np.radians(x/(3.6e6))) 
                        for x in np.power(10, ruwe_dist['log(sep)']) ]

# Adding: transform the reference data as well
ref_inputs = np.array([[float(row["log(sep)"]),float(row["DeltaG"])] 
                      for row in ref])
ref_inputs_lin = np.copy(ref_inputs)
ref_inputs_lin[:,0] = 10**ref_inputs[:,0]
ref_inputs_star = np.copy(ref_inputs)
ref_inputs_star[:,0] = dist_au * np.radians(ref_inputs_lin[:,0]/3.6e6)

#  2D interpolation functions for ruwe and sigma_ruwe
x_edges = np.unique(np.array(ruwe_dist['Sep(AU)']))
y_edges = np.unique(np.array(ruwe_dist['DeltaG']))
z = np.reshape(np.array(ruwe_dist['log(RUWE)']), 
               [len(y_edges), len(x_edges)])
z_sigma = np.reshape(np.array(ruwe_dist['sigma_log(RUWE)']), 
                     [len(y_edges), len(x_edges)])

f_ruwe2d = scipy.interpolate.interp2d(x_edges, y_edges, z)
f_sigma2d = scipy.interpolate.interp2d(x_edges, y_edges, z_sigma)

pred_log_ruwe2d = np.concatenate([f_ruwe2d(ref_inputs_star[i][0],
                                           ref_inputs_star[i][1]) 
                                           for i in range(nref)])
predicted_ruwe2d = 10**pred_log_ruwe2d
pred_sigma2d = np.concatenate([f_sigma2d(ref_inputs_star[i][0],
                                         ref_inputs_star[i][1]) 
                                         for i in range(nref)])


### Use the exact same input arrays, but interpolate with the RGI instead
f_ruwe = scipy.interpolate.RegularGridInterpolator(points=(x_edges, y_edges), 
                                                   values=z.T,
                                                   bounds_error=False,
                                                   fill_value=np.nan)
pred_log_ruwe = np.concatenate([f_ruwe(ref_inputs_star[i]) 
                                for i in range(nref)])



if __name__=="__main__":
    print(np.nanmax(pred_log_ruwe),np.nanmin(pred_log_ruwe))

    print(np.shape(pred_log_ruwe))
    print(np.shape(pred_log_ruwe2d))

    check_match = np.isclose(pred_log_ruwe,pred_log_ruwe2d,
                             atol=1e-7,equal_nan=True)

    print(len(np.where(check_match)[0]))

    for i in np.where(~check_match)[0]:
        print(ref_inputs_star[i],pred_log_ruwe2d[i],pred_log_ruwe[i])

def test_dist():
    assert pytest.approx(dist_au,rel=5e-4)==star_distance


# The different interpolators treat the lower separation bound differently
# interp2d gives a result and the RGI gives a nan
# So ignore them for testing purposes if everything else is fine
not_edge_cases = ref_inputs[:,0]>1.275
if __name__=="__main__":
    print("edge cases:",len(np.where(~not_edge_cases)[0]))
    print("mismatches that aren't edge cases",
          len(np.where(not_edge_cases & ~check_match)[0]))

def test_interp_methods():

    assert np.allclose(pred_log_ruwe[not_edge_cases],
                       pred_log_ruwe2d[not_edge_cases])

def test_comprehension():
    pred_log_ruwe_array = f_ruwe(ref_inputs_star[:,:2])

    assert np.allclose(pred_log_ruwe_array[not_edge_cases],
                       pred_log_ruwe2d[not_edge_cases])

# OK, so the issue is not the interpolation functions themselves, 
# we must have changed something about the ruwe code itself.
# The code above is now our benchmark for all future versions
# of ruwe.py

# Generate a set of companions for testing
comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
comps.generate()

ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
ruwe.read_dist()
f_ruwe_current, f_sigma_current = ruwe.interp_ruwe(dist_au)

pred_log_ruwe_current = f_ruwe(ref_inputs_star[:,:2])

def test_current_interp2d():
    assert np.allclose(pred_log_ruwe_current[not_edge_cases],
                       pred_log_ruwe2d[not_edge_cases])

def test_current_interp():
    assert np.allclose(pred_log_ruwe_current[not_edge_cases],
                       pred_log_ruwe[not_edge_cases])

