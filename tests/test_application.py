
import pytest

import pathlib, os, sys
import argparse

from molusc.application import Application
from molusc.companions import Companions
from molusc.utils import set_null_limits

limits = set_null_limits()
star_mass = 1.0
pd_mu = 5.03
pd_sig = 2.28
q_exp = 0.0

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
bfile = os.path.join(repo_path,"molusc/BinaryStarGUI_laptop.py")
cfile = os.path.join(repo_path,"example_contrast.txt")
rfile = os.path.join(repo_path,"example_rv.txt")

def cleanup(prefix):
    if os.path.exists(prefix+"_all.csv"):
        os.remove(prefix+"_all.csv")
    if os.path.exists(prefix+"_kept.csv"):
        os.remove(prefix+"_kept.csv")
    if os.path.exists(prefix+"_RVs.csv"):
        os.remove(prefix+"_RVs.csv")
    if os.path.exists(prefix+"_params_output.yml"):
        os.remove(prefix+"_params_output.yml")

test_fname = os.path.join(repo_path,"tests/test_write.hdf5")
# TODO: make a test companions file to use for read tests, 
# instead of writing it out every time
def test_write():
    """ Write out the test companions file for further use. """

    comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
    comps.generate()

    comps.write(test_fname)

def test_positional(monkeypatch):

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0], "cl", "mp_tst", 
                  "13h50m06.28s", "--",  "-40d50m08.9s", "100", "1.1"])
        print(sys.argv)
        app = Application(sys.argv)
        acheck = app.parse_input(app.input_args[1:])
        assert acheck


def test_optional(monkeypatch):
    # assert 0 == 0

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0], "cl", "--ao", cfile,  
                  "--filter", "K", "--comps", test_fname, 
                  "--",  "mp_tst", 
                  "13h50m06.28s",  "-40d50m08.9s", "100", 
                  "1.1"])
        app = Application(sys.argv)
        acheck = app.parse_input(app.input_args[1:])
        assert acheck

def test_infile(monkeypatch):
    # assert 0 == 0

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0], "cl", "--ao", cfile,  
                  "--filter", "K", "--comps", test_fname, 
                  "--",  "tests/mp_tst", 
                  "13h50m06.28s",  "-40d50m08.9s", "100", 
                  "1.1"])
        app = Application(sys.argv)
        acheck = app.parse_input(app.input_args[1:])
        assert app.companions_filename==test_fname

aprefix = os.path.join(repo_path,"tests/ao_tst")
def test_ao(monkeypatch):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0],"cl", "-a", "--ao", cfile,  
                  "--filter", "K", 
                  "--", aprefix, 
                  "13h50m06.28s", "-40d50m08.9s", "100", "1.1"])
        print(sys.argv)
        app = Application(sys.argv)
        app.start()
        assert os.path.exists(aprefix+"_all.csv")
    cleanup(aprefix)

rprefix = os.path.join(repo_path,"tests/rv_tst")
def test_rv(monkeypatch):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0],"cl", "-a", "--rv", rfile,  
                  "--resolution","50000", 
                  "--", rprefix, 
                  "13h50m06.28s", "-40d50m08.9s", "100", "1.1"])
        print(sys.argv)
        app = Application(sys.argv)
        app.start()
        assert os.path.exists(rprefix+"_RVs.csv")
    cleanup(rprefix)

