import os, sys, glob

import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import astropy.io.ascii as at
import h5py
import yaml

from utils import *

# rcParams['font.family'] = 'serif'
# rcParams['font.serif'] = ['Arial']
rcParams['font.weight'] = 'bold'
rcParams['axes.labelweight'] = 'bold'
rcParams['axes.linewidth'] = 1

# output_dir = os.get_env("DATA_PATH","./")
output_dir = "/Users/douglste/Google Drive/Shared drives/DouglasGroup/data"

def setup_axes():
    """
    Make a figure/axis object with appropriate limits
    """
    rcParams['xtick.labelsize'] = 'small'
    rcParams['ytick.labelsize'] = 'small'
    rcParams['axes.labelsize'] = 10

    fig, ax = plt.subplots(figsize=(5,3)) 
    ax.set_ylim(60, 700)
    ax.set_xlim(0.05,2e9)

    # Add a secondary axis, showing mass in solar masses
    secax = ax.secondary_yaxis('right', functions=(jup_mass_to_sol, sol_mass_to_jup))
    secax.set_ylabel(r'Mass ($M_\odot$)')
    # Add a secondary x axis, showing semi-major axis for a solar mass companion
    triax = ax.secondary_xaxis('top', functions=(period_to_a, a_to_period))
    triax.set_xlabel('a (AU)')
    triax.tick_params(axis='x', which='minor', bottom=False, top=False)
    # Correct axis scales, labels, and ticks
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Period (days)')
    ax.set_ylabel(r'Mass $(M_{Jup})$')
    plt.xticks()
    plt.yticks()

    return ax

def calc_limits(star,percentile):
    """
    Calculate the percentile limit for a selected star. 

    Currently, only one percentile value (a float or integer) can be provided
    """
    file_in = get_info(star)[0] #survivors
    star_mass = get_info(star)[5]
    # print(f"Star mass: {star_mass}")

    fig, ax = plt.subplots(figsize=(3.352242, 2.514181))  # this is sized to fit in one column of AAS journal format
    # Read in survivor data
    if file_in.endswith(".csv"):
        survivors = Table.read(file_in, format='ascii.csv')
    elif file_in.endswith(".h5") or file_in.endswith(".hdf5"):
        t = h5py.File(file_in,'r+')
        survivors = t["kept"]
    else:
        print("file type unknown:",file_in)
        return

    try:
        survivors['mass (M_Jup)'] = survivors['mass ratio']*star_mass*1047.35
    except:
        print("mass (M_Jup) array already exists!")

    try:
        survivors['logPeriod'] = np.log10(survivors['period(days)'])
    except:
        print("logPeriod array already exists!")

    # Need to split into logarithmic period bins 
    # and take the desired percentile of the masses in each bin
    bins = np.arange(0., 10., .5)
    # n, _ = np.histogram(np.log10(survivors['period(days)']), bins=bins)
    # Create a subchart for each bin. Find the 95 percentile mass in that bin, and the median period
    period = []
    perc = []
    # n = []
    for j in range(1, len(bins)):
        in_bin = (bins[j-1] <= survivors['logPeriod']) & (survivors['logPeriod'] <= bins[j])
        nbin = len(np.where(in_bin)[0])
        if nbin > 0:
            perc.append(np.percentile(survivors['mass (M_Jup)'][in_bin], percentile))
            period.append(np.median(survivors['period(days)'][in_bin]))
            # n.append(nbin) 

    return period, perc

def plot_limits(ax, period, perc, **plot_kwargs):
    """
    Plot the percentile limit for a selected star onto the provided axis
    """
    ax.plot(period, perc, **plot_kwargs)

def plot_detections():
    """
    For stars with AO detections, plot the raw projected separation as an arrow
    (I haven't calculated mass estimates yet)
    """
    pass

def choose_stars():
    """
    Run through the list of observed stars, and decide whether they should be 
    plotted as a limit or a detection
    """
    pass

def read_sample(data_dir=output_dir):
    """
    Open a target file and read in the list of observed stars
    """

    # Detection limits for everything observed
    # (doubles as the complete obs list)
    lim_columns = ["Name","MJD","Filter","Nobs","t_int",
                   "rho150","rho200","rho250","rho300","rho400","rho500",
                   "rho700","rho1000","rho1500","rho2000","PI"]
    keck_limfile = os.path.join(data_dir,"keck_psf_detections_praesepe/detlimtable_praesepe.txt")
    keck_alllim = at.read(keck_limfile,delimiter="&",data_start=1,
                         names=lim_columns)

    ksubidx = (keck_alllim["Filter"]=="Kc") | (keck_alllim["Filter"]=="Kp")    
    keck_douglas = keck_alllim[(keck_alllim["PI"]=="Douglas") & ksubidx]
    unames, cts = np.unique(keck_douglas["Name"], return_counts=True)

    # If there are any duplicated rows, choose the one with the deepest contrast
    if np.any(cts>1):
        for ui in np.where(cts>1)[0]:
            bad_name = unames[ui]
            bad_loc = np.where(keck_douglas["Name"]==bad_name)[0]
            check_row = []
            for colname in lim_columns:
                if ("rho" in colname)==False:
                    continue
                check_row = np.argmax(keck_douglas[colname][bad_loc])
            rvals, rcts = np.unique(check_row, return_counts = True)
            best_row = rvals[np.argmax(rcts)]
            keep_row = bad_loc[best_row]
            del_row = np.delete(bad_loc,bad_loc==keep_row)

            keck_douglas.remove_rows(del_row)

    unames, cts = np.unique(keck_douglas["Name"], return_counts=True)
    if np.any(cts>1):
        print("DUPLICATE REMOVAL FAILED")


    print(len(keck_douglas),"observations")

    new_keck = []
    keck_cand = []

    # Adam's direct imaging detections
    det_columns = ["Name","MJD","Filter","Nobs","rho","PA","Deltam","PI"]
    keck_detfile = os.path.join(data_dir,"keck_psf_detections_praesepe/dettable_praesepe.txt")
    keck_alldet = at.read(keck_detfile,delimiter="&",data_start=1,
                         names=det_columns)
    subnames,subidx = np.unique(keck_alldet["Name"],return_index=True)
    keck_alldet = keck_alldet[subidx]
    keck_alldet["cat_idx"] = np.ones(len(keck_alldet),int)*-99
    # print(keck_alldet.dtype)

    for i,newname in enumerate(keck_alldet["Name"]):
        loc = np.where(keck_douglas["Name"]==newname)[0]

        if len(loc)==1:
            # print(newname,loc,pdat["BINARY"][loc])
            keck_alldet["cat_idx"][i] = loc[0]
            sep = np.float32(keck_alldet["rho"][i].split("$")[0])
            if (sep<2000):
                new_keck.append(newname.replace(" ","_"))
                keck_alldet["cat_idx"][i] = loc[0]
            elif (sep>=2000):
                keck_cand.append(newname.replace(" ","_"))
                keck_alldet["cat_idx"][i] = loc[0]
        elif len(loc)>1:
            print("Duplicate!",newname,loc)
            print(keck_douglas[loc])
        else:
            continue
            # print("Not found",newname,loc)

    # Adam's PSF-fitting results
    mult_columns = ["Name","MJD","Filter","Nobs","rho","PA","Deltam","PI"]
    keck_multfile = os.path.join(data_dir,"keck_psf_detections_praesepe/multitable_praesepe.txt")
    keck_allmult = at.read(keck_multfile,delimiter="&",data_start=1,
                         names=mult_columns)
    ksubidx = (keck_allmult["Filter"]=="Kc") | (keck_allmult["Filter"]=="Kp")
    keck_allmult = keck_allmult[ksubidx]

    subnames,subidx = np.unique(keck_allmult["Name"],return_index=True)
    keck_allmult = keck_allmult[subidx]
    # keck_allmult.dtype
    keck_allmult["cat_idx"] = np.ones(len(keck_allmult),int)*-99
    for i,newname in enumerate(keck_allmult["Name"]):

        loc = np.where(keck_douglas["Name"]==newname)[0]
        if len(loc)==1:
            new_keck.append(newname.replace(" ","_"))
            keck_allmult["cat_idx"][i] = loc[0]
        elif len(loc)>1:
            print("Duplicate!",newname,loc)
            print(keck_douglas[loc])
        else:
            continue
            # print("Not found",newname,loc)

    # Aaron's masking results
    # (includes some duplicates of Adam's PSF-fitting detections)
    mask_file = os.path.join(data_dir,"keck_masking_detections_praesepe/masking_detections.csv")
    keck_mask = at.read(mask_file)
    keck_mask["cat_idx"] = np.ones(len(keck_mask),int)*-99
    # print(keck_mask.dtype)

    for i,newname in enumerate(keck_mask["Name"][(keck_mask["Sep"]!="SINGLE") & (keck_mask["Sep"]!="no cals")]):
        loc = np.where(keck_douglas["Name"]==newname)[0]

        if (len(loc)==1) and ((newname in np.append(new_keck,keck_cand))==False):
            keck_mask["cat_idx"][i] = loc[0]
            # print(newname,loc,pdat["BINARY"][loc])
            new_keck.append(newname.replace(" ","_"))
        elif (len(loc)>1):
            print("Duplicate!",newname,loc)
            print(keck_douglas[loc])
        else:
            continue
            # print("Not found",newname,loc)


    new_keck = np.array(new_keck)
    print(new_keck)
    keck_cand = np.array(keck_cand)
    print(keck_cand)

    all_new = np.unique(np.append(new_keck,keck_cand))
    # Remove Andrew's targets
    all_new = np.delete(all_new,all_new=="PM_I08131-1355")
    all_new = np.delete(all_new,all_new=="PM_I10367+1521")
    all_new = np.delete(all_new,all_new=="HIP14807")
    all_new = np.delete(all_new,all_new=="BD+23__635")
    all_new = np.delete(all_new,all_new=="03562+5939")

    print(len(all_new))

    return all_new, keck_douglas, keck_alldet, keck_allmult, keck_mask

if __name__=="__main__":
    all_new, keck_douglas, keck_alldet, keck_allmult, keck_mask = read_sample()

    ax = setup_axes()
    
    for name in all_new:
        try:
            period, perc = calc_limits(name,90)
            plot_limits(ax, period, perc)
        except:
            print(name,"results not found")
            continue

    plt.show()