import pathlib, os
import warnings

import h5py
import pytest
import numpy as np
import astropy.io.ascii as at
import scipy.stats as stats

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
from molusc.gaia import Gaia
from molusc.companions import Companions
from molusc.utils import set_null_limits

warnings.simplefilter('error', category=RuntimeWarning)


limits = set_null_limits()
pd_mu = 5.03
pd_sig = 2.28
q_exp = 0.0

# JS178//Gaia DR3 659777074828577024, which has a contaminant in Gaia
# Nearest neighbor is Gaia DR3 659777070533194112
# 2.775 arcsec away by Vizier’s calculations
star_ra = '08h37m24.18s '
star_dec = '+19d25m1.2s'
star_age = 0.8
star_mass = 0.597
plx = 5.387 #mas
star_ruwe = 1.0249 #Rizzuto+ 2020

dist_pc = 1000/plx 
dist_au = dist_pc * 206265

star_distance = 1 / (plx) * 2.063e+8 # AU - Mackenna's calculation

nn_sep_as = 2.775
nn_sep_deg = nn_sep_as/3600
nn_sep_rad = np.radians(nn_sep_deg)
nn_sep_au = dist_au * np.tan(nn_sep_rad)

gfile = os.path.join(repo_path, "reference_data/gaia_contrast.txt")
# ref = at.read(ref_file)
# ref_inputs = np.array([[float(row["log(sep)"]),float(row["DeltaG"])] for row in ref])

# Generate a set of companions for testing
comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
comps.generate()

gaia_limit=18

def test_gaia_init():
    gg = Gaia(gfile, comps, star_mass, star_age, star_ra, 
              star_dec, "K", gaia=True)
    assert ((star_ra==gg.star_ra) and (star_dec==gg.star_dec)
            and (star_mass==gg.star_mass) and (gg.a_type=="gaia") )

def test_gaia_parallax():
    gg = Gaia(gfile, comps, star_mass, star_age, star_ra, 
              star_dec, "K", gaia=True)
    dist_fail = gg.get_distance(plx)
    assert ((dist_fail==0) and 
            (pytest.approx(dist_au,rel=1e-4)==gg.star_distance))

def test_gaia_nn_distance():
    # Just making sure it runs
    gg = Gaia(gfile, comps, star_mass, star_age, star_ra, 
              star_dec, "K", gaia=True)
    dist_fail = gg.get_distance(plx,get_neighbor=True)
    assert (((dist_fail==0) or (dist_fail==-52))and 
            (pytest.approx(gg.nearest_neighbor_dist,rel=1e-4)==(nn_sep_au)))

def test_gaia_read():
    gg = Gaia(gfile, comps, star_mass, star_age, star_ra, 
              star_dec, "K", gaia=True)
    read_fail = gg.read_contrast()
    assert (read_fail==0) and (gg.a_type!='hard limit')

def test_gaia_run():
    # Just making sure it runs
    gg = Gaia(gfile, comps, star_mass, star_age, star_ra, 
              star_dec, "K", gaia=True)
    dist_fail = gg.get_distance(plx,get_neighbor=True)
    read_fail = gg.read_contrast()
    res = gg.analyze_gaia(gaia_limit)
    assert res is not None

def test_neighbors():
    gg = Gaia(gfile, comps, star_mass, star_age, star_ra, 
              star_dec, "K", gaia=True)
    dist_fail = gg.get_distance(plx,get_neighbor=True)
    read_fail = gg.read_contrast()
    res = gg.analyze_gaia(gaia_limit,nearest_neighbor=True)
    assert (np.all(gg.pro_sep[res==False]<nn_sep_au))

