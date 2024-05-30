
import pytest

import pathlib, os, sys
import argparse
import h5py

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
    if os.path.exists(prefix+"_all.h5"):
        os.remove(prefix+"_all.h5")
    if os.path.exists(prefix+"_kept.h5"):
        os.remove(prefix+"_kept.h5")
    if os.path.exists(prefix+".h5"):
        os.remove(prefix+".h5")
    if os.path.exists(prefix+"_RVs.h5"):
        os.remove(prefix+"_RVs.h5")
    if os.path.exists(prefix+"_params_output.yml"):
        os.remove(prefix+"_params_output.yml")

# If the number of companions is above a certain number, then 
# the file written out should be an hdf5 file, not csv
aprefix = os.path.join(repo_path,"tests/ao_tst")
def test_h5write(monkeypatch):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0],"cl", "-a", "--ao", cfile,  
                  "--filter", "K", 
                  "--", aprefix, 
                  "13h50m06.28s", "-40d50m08.9s", "50_001", "1.1"])
        print(sys.argv)
        app = Application(sys.argv)
        app.start()
        assert os.path.exists(aprefix+".h5")

def test_h5post(monkeypatch):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0],"cl", "--ao", cfile,  
                  "--filter", "K", 
                  "--", aprefix, 
                  "13h50m06.28s", "-40d50m08.9s", "50_001", "1.1"])
        print(sys.argv)
        app = Application(sys.argv)
        app.start()
        with h5py.File(aprefix+".h5","r") as f:
            assert "kept" in f.keys()

def test_h5noprior(monkeypatch):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0],"cl", "--ao", cfile,  
                  "--filter", "K", 
                  "--", aprefix, 
                  "13h50m06.28s", "-40d50m08.9s", "50_001", "1.1"])
        print(sys.argv)
        app = Application(sys.argv)
        app.start()
        with h5py.File(aprefix+".h5","r") as f:
            assert ("all" in f.keys())==False

def test_h5prior(monkeypatch):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0],"cl", "-a", "--ao", cfile,  
                  "--filter", "K", 
                  "--", aprefix, 
                  "13h50m06.28s", "-40d50m08.9s", "50_001", "1.1"])
        print(sys.argv)
        app = Application(sys.argv)
        app.start()

        with h5py.File(aprefix+".h5","r") as f:
            assert "all" in f.keys()
    cleanup(aprefix)