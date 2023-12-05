import os, sys, pathlib

from molusc.companions import Companions

default_cache = os.expanduser("~/")
MOLUSC_CACHE_PATH = os.getenv("MOLUSC_CACHE_PATH",default_cache)

if __name__=="__main__":
    # Using the default parameters from application.py for now
    num_generated = 50#_000_000
    limits = [None]*21
    pd_mu = 5.03
    pd_sig = 2.28
    q_exp = 0.0

    comps = Companions(num_generated, limits, pd_mu, pd_sig, q_exp)
    comps.generate()


    comps.write(os.path.join(MOLUSC_CACHE_PATH,"MOLUSC_prior_test.hdf5"))
