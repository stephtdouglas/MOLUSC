import pathlib, os

import h5py
import pytest
import numpy as np

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
from molusc.companions import Companions

limits = [None]*21
star_mass = 1.0
pd_mu = 5.03
pd_sig = 2.28
q_exp = 0.0

test_fname = os.path.join(repo_path,"tests/test_write.hdf5")

def test_generation():
    """Ensure companion generation executes cleanly."""

    comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
    comps.generate()

def test_write():
    """ Ensure the write function executes cleanly. """

    comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
    comps.generate()

    comps.write(test_fname)

def test_write_values():
    """ Ensure the write function executes cleanly. """

    comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
    comps.generate()

    comps.write(test_fname)
    with h5py.File(test_fname,"r") as f:

        # If all the period values are still zeros, that's a problem
        are_periods_zero = f["companions"][0] == 0
        assert np.all(are_periods_zero)==False

def test_bad_extension():
    """ 
    Ensure the write function raises an error with the wrong 
    file extension 
    """

    comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)

    test_fname2 = os.path.join(repo_path,"tests/test_write.csv")

    with pytest.raises(ValueError):
        comps.write(test_fname2)

def test_read():
    """ Ensure the write function executes cleanly. """

    comps = Companions.read(test_fname)