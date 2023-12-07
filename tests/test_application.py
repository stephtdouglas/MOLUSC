
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

test_fname = os.path.join(repo_path,"tests/test_write.hdf5")
# TODO: make a test companions file to use for read tests, 
# instead of writing it out every time
def test_write():
    """ Write out the test companions file for further use. """

    comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
    comps.generate()

    comps.write(test_fname)

def test_positional(monkeypatch):
    # assert 0 == 0

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0], "cl", "mp_test", 
                  "13h50m06.28s", "--",  "-40d50m08.9s", "100", "1.3"])
        print(sys.argv)
        app = Application(sys.argv)
        acheck = app.parse_input(app.input_args[1:])
        assert acheck
        # assert False


def test_optional(monkeypatch):
    # assert 0 == 0

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0], "cl", "--ao", cfile,  
                  "--filter", "K", "--comps", test_fname, 
                  "--",  "mp_test", 
                  "13h50m06.28s",  "-40d50m08.9s", "100", 
                  "1.3"])
        app = Application(sys.argv)
        acheck = app.parse_input(app.input_args[1:])
        assert acheck

def test_infile(monkeypatch):
    # assert 0 == 0

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0], "cl", "--ao", cfile,  
                  "--filter", "K", "--comps", test_fname, 
                  "--",  "mp_test", 
                  "13h50m06.28s",  "-40d50m08.9s", "100", 
                  "1.3"])
        app = Application(sys.argv)
        acheck = app.parse_input(app.input_args[1:])
        assert app.companions_filename!=""

