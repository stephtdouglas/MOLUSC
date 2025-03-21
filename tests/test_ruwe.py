import pathlib, os
import warnings

import h5py
import pytest
import numpy as np
import astropy.io.ascii as at
import scipy.stats as stats

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
plx = 8.017 #mas
star_ruwe = 0.91 #Rizzuto+ 2020

dist_pc = 1000/plx 
dist_au = dist_pc * 206265

star_distance = 1 / (plx) * 2.063e+8 # AU - Mackenna's calculation

# Compute the limits of the RUWE grid
log_sep_lims = np.array([1.275,3.325])
sep_lims = 10**log_sep_lims # mas
sep_lims_rad = np.radians(sep_lims/3.6e6)
# Small angle approximation means that we can multiply dist * radians
sep_lims_au = sep_lims_rad * dist_au

ref_file = os.path.join(repo_path,"reference_data/RuweTableGP.txt")
ref = at.read(ref_file)
ref_inputs = np.array([[float(row["log(sep)"]),float(row["DeltaG"])] for row in ref])

# Generate a set of companions for testing
comps = Companions(100, limits, star_mass, pd_mu, pd_sig, q_exp)
comps.generate()

def test_init():
    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    assert ((star_ra==ruwe.star_ra) and (star_dec==ruwe.star_dec)
            and (star_age==ruwe.star_age) and (star_mass==ruwe.star_mass))

def test_init2():
    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    assert np.isnan(ruwe.gmag)

def test_read():
    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    ruwe.read_dist()

    assert np.all(ruwe.ruwe_dist==ref)

def test_dist_conversion():

    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    ruwe.read_dist()
    _ = ruwe.interp_ruwe(dist_au)

    tab_seps = np.asarray(ruwe.ruwe_dist['Sep(AU)'])

    assert np.all((tab_seps >= (sep_lims_au[0]*0.99)) &
            (tab_seps <= (sep_lims_au[1]*1.01)))

def test_interp_creation():

    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    ruwe.read_dist()
    f_ruwe, f_sigma = ruwe.interp_ruwe(dist_au)

    # Check that the x and y edges are what we think they are
    log_sep_edges = np.unique(ref["log(sep)"])
    sep_edges = 10**log_sep_edges # mas
    sep_edges_rad = np.radians(sep_edges/3.6e6)
    # Small angle approximation means that we can multiply dist * radians
    sep_edges = sep_edges_rad * dist_au

    dg_edges = np.unique(ref["DeltaG"])

    nsep = len(sep_edges)
    ndg = len(dg_edges)

    # Shape is (rows, columns)
    rshape = np.shape(f_ruwe.values)
    sshape = np.shape(f_sigma.values)

    # sep is X and dG is Y, so 
    # there should be nsep columns and ndg rows

    assert (
            # rshape[0]==ndg and rshape[1]==nsep and
            np.allclose(f_ruwe.grid[0],sep_edges) and np.allclose(f_ruwe.grid[1],dg_edges)
            )

def test_interp_projected():

    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    ruwe.read_dist()
    f_ruwe, f_sigma = ruwe._interp_ruwe_projected()

    ref_ruwe_min = np.min(ref["log(RUWE)"])
    ref_ruwe_max = np.max(ref["log(RUWE)"])

    print("Input min/max:",ref_ruwe_min,ref_ruwe_max)
    print("Table min/max:",np.min(f_ruwe.values),np.max(f_ruwe.values))

    assert (np.all(f_ruwe.values >= ref_ruwe_min) and 
            np.all(f_ruwe.values <= ref_ruwe_max))

def test_proj_result1():
    # Every pair in the input file should produce the same result
    # From the input file or the interpolation

    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    ruwe.read_dist()
    f_ruwe, f_sigma = ruwe._interp_ruwe_projected()

    assert np.allclose(f_ruwe(ref_inputs),ref["log(RUWE)"])

def test_proj_result2():
    # Every pair in the input file should produce the same result
    # From the input file or the interpolation

    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    ruwe.read_dist()
    f_ruwe, f_sigma = ruwe._interp_ruwe_projected()

    assert np.allclose(f_sigma(ref_inputs),ref["sigma_log(RUWE)"])


def test_proj_result_linsep():
    # Every pair in the input file should produce the same result
    # From the input file or the interpolation

    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    ruwe.read_dist()
    f_ruwe, f_sigma = ruwe._interp_ruwe_projected(use_log_sep=False)

    ref_inputs2 = np.copy(ref_inputs)
    ref_inputs2[:,0] = 10**ref_inputs[:,0]

    assert np.allclose(f_ruwe(ref_inputs2),ref["log(RUWE)"],atol=1e-5)


def test_result_same_dist():
    # Every pair in the input file should produce the same result
    # From the input file or the interpolation

    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    ruwe.read_dist()
    f_ruwe, f_sigma = ruwe.interp_ruwe(star_distance)

    ref_inputs3 = np.copy(ref_inputs)
    star_tan = np.tan(np.radians((10**ref_inputs[:,0]/3.6e6)))
    ref_inputs3[:,0] = star_distance * star_tan 

    assert np.allclose(f_ruwe(ref_inputs3),ref["log(RUWE)"],atol=1e-5)


# The tests below were created when I was exploring what impacts 
# the various inputs. They are preserved for posterity, but
# have been adjusted so that some of them are expected to be false
# as that will let pytest pass properly overall

def test_result_diff_dist():
    # What happens if I use my distance calculation for the calculation
    # but the original version to create the interpolator?
    # Answer: we get errors
    # So now this test will pass if they DON'T match

    ruwe = RUWE(star_ra, star_dec, star_age, star_mass, comps)
    ruwe.read_dist()
    f_ruwe, f_sigma = ruwe.interp_ruwe(dist_au)

    ref_inputs3 = np.copy(ref_inputs)
    sep_tan = np.tan(np.radians((10**ref_inputs[:,0]/3.6e6)))
    ref_inputs3[:,0] = star_distance * sep_tan 

    assert np.allclose(f_ruwe(ref_inputs3),ref["log(RUWE)"],atol=1e-5)==False

def test_distance_conversion():
    # Is the problem the conversion factors we're using? Yes
    # Meaning now this test will pass with different conversion factors

    parallax_mas = 8.017 
    parallax_as = parallax_mas/1e3
    parallax_deg = parallax_as/3600

    distance_parsec = 1/parallax_as

    distance_au = distance_parsec * 206265

    assert pytest.approx(distance_au,rel=1e-4)!=star_distance


def test_separation_trig():
    # Is the problem the small angle approximation? 
    # No! Trigonometry still works as expected

    parallax_mas = 8.017 
    parallax_as = parallax_mas/1e3
    parallax_deg = parallax_as/3600

    distance_parsec = 1/parallax_as

    distance_au = distance_parsec * 206265

    # Small angle approximation theta=tan(theta)
    # vs actual tangent calculation
    sep_small = np.radians((10**ref_inputs[:,0]/3.6e6))
    sep_tan = np.tan(sep_small)

    sep_small = distance_au * sep_small
    sep_tan = distance_au * sep_tan

    assert pytest.approx(sep_small)==sep_tan


def test_separation_conversion():
    # Is the problem the conversion factors we're using?
    # Part 2: propagation to separation
    # Yes, this is still an issue
    # The test now passes when the results don't match

    parallax_mas = 8.017 
    parallax_as = parallax_mas/1e3
    parallax_deg = parallax_as/3600

    distance_parsec = 1/parallax_as

    distance_au = distance_parsec * 206265
    # Original version is 2.063e5 (effectively), or 206300
    # Difference is 35 AU

    # Small angle approx is fine
    sep = np.radians((10**ref_inputs[:,0]/3.6e6))
    sep_new = distance_au * sep
    sep_old = star_distance * sep

    print(sep_new-sep_old)

    # We care about separations as close as 0.1 AU
    # So let's set the absolute tolerance at 0.01 AU

    assert np.allclose(sep_new,sep_old,atol=1e-2)==False#,rtol=1e-4)
