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
# plx = 8.017 #mas

# Generate a set of companions for testing
comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
comps.generate()


ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
ruwe.read_dist()
f_ruwe, f_sigma = ruwe._interp_ruwe_projected(use_log_sep=False)

# sep (mas), delta-g, expected minimum logRUWE
# -99 is for when I expect a NaN
pytestmark = pytest.mark.parametrize('set_min',
                                    [[300,0.5,1.5],
                                    [9,1.5,0.7],
                                    [100,3.5,0.3],
                                    [3,6,0.01],
                                    [3,-2,-99],
                                    # These specify a maximum, so we expect the tests to fail
                                    pytest.param([60,0.5,1.6],marks=[pytest.mark.xfail(strict=True)]),
                                    pytest.param([2000,1,0.5],marks=[pytest.mark.xfail(strict=True)]),
                                    # These should result in a nan from the interpolation function
                                    pytest.param([10,6,0.1]),
                                    pytest.param([1,8,-99]),
                                    ])

def test_interp_result_min(set_min):

    res = f_ruwe(set_min[0:2])
    print(set_min, res)
    assert np.isnan(res) or (set_min[2] <= res) 

# def calc_rejection_prob(f_ruwe, f_sigma, sg):
#     res = 10**f_ruwe(sg)
#     sig = 10**f_sigma(sg)
#     rejection_prob = stats.halfnorm.cdf(res,loc=star_ruwe, scale=10**sig)
#     return rejection_prob   
