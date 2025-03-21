import pathlib, os
import warnings

import h5py
import pytest
import numpy as np
import astropy.io.ascii as at
import scipy.stats as stats

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
from molusc.ao import AO
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

star_distance = 1 / (plx) * 2.063e+8 # AU - Mackenna's calculation

ex_file = os.path.join(repo_path,"example_contrast.txt")
# ref = at.read(ref_file)
# ref_inputs = np.array([[float(row["log(sep)"]),float(row["DeltaG"])] for row in ref])

# Generate a set of companions for testing
comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
comps.generate()

def test_init():
    ao = AO(ex_file, comps, star_mass, star_age, star_ra, 
              star_dec, "K", gaia=False)
    assert ((star_ra==ao.star_ra) and (star_dec==ao.star_dec)
            and (star_mass==ao.star_mass)
            and (ao.ao_filename==ex_file) and (ao.a_type=="") )

def test_parallax():
    ao = AO(ex_file, comps, star_mass, star_age, star_ra, 
              star_dec, "K", gaia=False)
    dist_fail = ao.get_distance(plx)
    assert ((dist_fail==0) and 
            (pytest.approx(dist_au,rel=1e-4)==ao.star_distance))

