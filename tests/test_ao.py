import pathlib, os

import h5py
import pytest
import numpy as np

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
from molusc.companions import Companions
from molusc.utils import set_null_limits
from molusc.ao import AO, AO_detection

limits = set_null_limits()
star_mass = 1.0
pd_mu = 5.03
pd_sig = 2.28
q_exp = 0.0
star_age = 0.8
star_ra = "13h50m06.28s" 
star_dec = "-40d50m08.9s"
test_filter = "K"

test_fname = os.path.join(repo_path,"tests/test_write.hdf5")
comps = Companions.read(test_fname)

test_aofile = os.path.join(repo_path,"example_contrast.txt")


def test_basic():
    assert True


def test_init():
    ao = AO(test_aofile,comps,star_mass,star_age,test_filter,gaia=False)

    assert ao.star_mass==star_mass


def test_dist():
    ao = AO(test_aofile,comps,star_mass,star_age,test_filter,gaia=False)

    ao.get_distance(star_ra, star_dec)


def test_cont():
    ao = AO(test_aofile,comps,star_mass,star_age,test_filter,gaia=False)

    ao.get_distance(star_ra, star_dec)
    ao.read_contrast()

def test_analyze():
    ao = AO(test_aofile,comps,star_mass,star_age,test_filter,gaia=False)

    ao.get_distance(star_ra, star_dec)
    ao.read_contrast()
    ao.analyze()

# EPIC 211998192 & 58456.59 &   Kp     &   2 &  3925.9 $\pm$ 7.4 &  148.770 $\pm$ 0.108 &  
# 6.314 $\pm$0.142 &              Douglas \\ % M N2.20181204.51037

def test_init_det():
    ao = AO_detection(comps,star_mass,star_age,"K",
                      3925.9, 7.4, #rho in mas
                      148.770, 0.108, # pa in deg
                      6.314, 0.142, # deltaM in mag
                      obs_date=58456.59)

    assert ao.star_mass==star_mass
