import os, sys, pathlib
import numpy as np

from molusc.companions import Companions
from molusc.utils import set_null_limits

default_cache = os.path.expanduser("~/molusc_cache/")
MOLUSC_CACHE_PATH = os.getenv("MOLUSC_CACHE_PATH",default_cache)

if __name__=="__main__":
    # Using the default parameters from application.py for now
    num_generated = 50_000_000
    limits = set_null_limits()
    star_mass = 0.5
    pd_mu = 5.03
    pd_sig = 2.28
    q_exp = 0.0

    limits["v0"]["mu"] = 35. #km/s
    limits["v0"]["sigma"] = 10. #km/s
    # limits["v0"]["fixed"]=35. #km/s
    # limits["v0"]["mu"] = None
    # limits["v0"]["sigma"] = None
    print(limits["v0"])

    limits["P"]["shape"] = "logflat"

    np.random.seed(424242)
    comps = Companions(num_generated, limits, star_mass, pd_mu, pd_sig, q_exp)
    comps.generate()

    comps.write(os.path.join(MOLUSC_CACHE_PATH,"MOLUSC_prior_v0_Pflat_50M.hdf5"))
