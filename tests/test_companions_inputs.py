import pathlib, os

import h5py
import pytest
import numpy as np

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
from molusc.companions import Companions
from molusc.utils import set_null_limits

pytestmark = pytest.mark.parametrize('P_fixed,a_fixed,q_fixed',
                                      [(36.525,None,None),
                                       (None,0.262074,None),
                                       (None,None,0.8),
                                       (36.525,0.262074,None),
                                       (None,0.262074,0.8),
                                       (36.525,None,0.8)
                                       ])

test_fname = os.path.join(repo_path,"tests/test_write.hdf5")

def test_fixed_paq2(P_fixed, q_fixed, a_fixed):
    """ make sure the code runs cleanly with various pairs of fixed params
    """

    bad_limits = set_null_limits()
    bad_limits["P"]["fixed"] = P_fixed
    bad_limits["q"]["fixed"] = q_fixed
    bad_limits["a"]["fixed"] = a_fixed
    star_mass = 1.0
    pd_mu = 5.03
    pd_sig = 2.28
    q_exp = 0.0

    comps = Companions(100, bad_limits, star_mass, pd_mu, pd_sig, q_exp)

    comps.generate()

    check_vals = True

    if P_fixed is not None:
        this_check = np.all(comps.P == pytest.approx(P_fixed))
        if this_check==False:
            print("P fixed failed")
        check_vals = check_vals & this_check

    if a_fixed is not None:
        this_check = np.all(comps.a == pytest.approx(a_fixed))
        if this_check==False:
            print("a fixed failed")
        check_vals = check_vals & this_check

    if q_fixed is not None:
        this_check = np.all(comps.mass_ratio == pytest.approx(q_fixed))
        if this_check==False:
            print("q fixed failed")
        check_vals = check_vals & this_check

    assert check_vals

def test_fixed_pairs_again(P_fixed, q_fixed, a_fixed):
    """ make sure the code runs cleanly with various pairs of fixed params
    """

    bad_limits = set_null_limits()
    bad_limits["P"]["fixed"] = P_fixed
    bad_limits["q"]["fixed"] = q_fixed
    bad_limits["a"]["fixed"] = a_fixed
    star_mass = 1.0
    pd_mu = 5.03
    pd_sig = 2.28
    q_exp = 0.0

    comps = Companions(100, bad_limits, star_mass, pd_mu, pd_sig, q_exp)

    comps.generate()

    if P_fixed is not None and q_fixed is not None:
        assert comps.a == pytest.approx(0.262074,rel=1e-3)

    elif a_fixed is not None and q_fixed is not None:
        assert comps.P == pytest.approx(36.525,rel=1e-3)

    elif P_fixed is not None and a_fixed is not None:
        assert comps.mass_ratio == pytest.approx(0.8, rel=1e-2)

    else:
        assert True