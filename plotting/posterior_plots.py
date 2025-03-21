#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, sys, glob, json

import matplotlib.pyplot as plt
import numpy as np
import h5py
import tables as tb
#import json,yaml

from sample_limits import setup_axes, calc_limits

data_path = os.getenv("DATA_PATH","/Users/douglste/data2/molusc_outputs")
output_dir = os.path.join(data_path,"molusc_outputs/")
print(output_dir)

#output_dir = "/scratch/douglste_JS117/"

def calc_limits2(res, mass_ratio, star_mass, period, colnames, percentile=95):
    
    # Need to split into logarithmic period bins 
    # and take the desired percentile of the masses in each bin
    bins = np.arange(0., 10., .5)
    # Create a subchart for each bin. Find the desired percentile mass in that bin, and the median period

    # Select which stars to include. We want to keep things that were not rejected.
    if isinstance(colnames,str):
        select = np.where(res["all"][colnames][:]==0)[0]
    else:
        select0 = np.ones(len(mass_ratio),bool)
        for colname in colnames:
            #print(colname,len(np.where(res["all"][colname])[0]))
            select0 = select0 & (res["all"][colname][:]==0)
        select = np.where(select0)[0]

    mass_msun = mass_ratio[select]*star_mass
    logperiod = np.log10(period[select])

    bin_idx = np.digitize(logperiod,bins=bins,right=False)
    
    nbins = len(bins)-1
    per_out, perc = np.zeros(nbins), np.zeros(nbins)
    
    for i in range(nbins):
        in_bin = np.where(bin_idx==i)[0]
        nbin = len(in_bin)
#        print(i,nbin)
        if nbin > 0:
            perc[i] = np.percentile(mass_msun[in_bin], percentile)
            per_out[i] = np.median(period[select][in_bin])
#        print("done",i)

    return per_out, perc


def calc_limits_q(res, mass_ratio, period, colnames, percentile=95):
    
    # Need to split into logarithmic period bins 
    # and take the desired percentile of the masses in each bin
    bins = np.arange(0., 10., .5)
    # Create a subchart for each bin. Find the desired percentile mass in that bin, and the median period

    # Select which stars to include. We want to keep things that were not rejected.
    if isinstance(colnames,str):
        select = np.where(res["all"][colnames][:]==0)[0]
    else:
        select0 = np.ones(len(mass_ratio),bool)
        for colname in colnames:
            #print(colname,len(np.where(res["all"][colname])[0]))
            select0 = select0 & (res["all"][colname][:]==0)
        select = np.where(select0)[0]

    mass_rat = mass_ratio[select]
    logperiod = np.log10(period[select])

    bin_idx = np.digitize(logperiod,bins=bins,right=False)
    
    nbins = len(bins)-1
    per_out, perc = np.zeros(nbins), np.zeros(nbins)
    
    for i in range(nbins):
        in_bin = np.where(bin_idx==i)[0]
        nbin = len(in_bin)
#        print(i,nbin)
        if nbin > 0:
            perc[i] = np.percentile(mass_rat[in_bin], percentile)
            per_out[i] = np.median(period[select][in_bin])
#        print("done",i)

    return per_out, perc


def plot_one(filename):
#    filename = os.path.join(output_dir,h5file)
    #print(os.path.exists(filename))
    res0 = tb.open_file(filename, drive="H5FD_CORE")
    res = res0.root
    #print(res)
    
    
    with h5py.File(filename,"r") as hf:
        meta = json.loads(hf["meta"]["yaml"].asstr()[()])
    star_mass = meta["star"]["mass"]
    
    #select = np.where((res["all"]["AO Rejected 1"][:]==0))[0]
    ax = setup_axes()
    for percentile in [90]:
        per_ao, perc_ao = calc_limits2(res,res["all"]['mass ratio'][:],
                                       star_mass,
                                       res["all"]['period(days)'][:],
                                       ["AO Rejected 1"],percentile=percentile)
        ax.plot(per_ao,np.asarray(perc_ao)*1047.35,'o')
        per_ao, perc_ao = calc_limits2(res,res["all"]['mass ratio'][:],
                                       star_mass,
                                       res["all"]['period(days)'][:],
                                       ["AO Rejected 1","Gaia Rejected"],percentile=percentile)
        ax.plot(per_ao,np.asarray(perc_ao)*1047.35,'o')
        per_ao, perc_ao = calc_limits2(res,res["all"]['mass ratio'][:],
                                       star_mass,
                                       res["all"]['period(days)'][:],
                                       ["AO Rejected 1","Gaia Rejected","RV Rejected"],percentile=percentile)
        ax.plot(per_ao,np.asarray(perc_ao)*1047.35,'o')
    plt.savefig(filename.replace(".h5",".png"))
    
    res0.close()


def plot_all(file_list,output_filename,percentile=90):

    ax = setup_axes("q")
    col_combos = [["AO Rejected 1"],["AO Rejected 1","Gaia Rejected"],
                  ["AO Rejected 1","Gaia Rejected","RV Rejected"]]
    colors = ["C0","C1","C2"]
    labels = ["AO only","AO+Gaia","AO+Gaia+RV"]
    for i in range(len(col_combos)):
        plt.plot([],[],"-",color=colors[i],alpha=0.75,label=labels[i])
    plt.legend(loc=4)

    for filename in file_list:
        #print(os.path.exists(filename))
        res0 = tb.open_file(filename, drive="H5FD_CORE")
        res = res0.root
        
        
        with h5py.File(filename,"r") as hf:
            meta = json.loads(hf["meta"]["yaml"].asstr()[()])
        star_mass = meta["star"]["mass"]

        for i,cc in enumerate(col_combos):
            per_ao, perc_ao = calc_limits_q(res,res["all"]['mass ratio'][:],
                                           #star_mass,
                                           res["all"]['period(days)'][:],
                                           cc,percentile=percentile)
            ax.plot(per_ao,np.asarray(perc_ao),#*1047.35,
                    '-',color=colors[i],
                    alpha=0.5,zorder=i*50)

        res0.close()

    plt.savefig(output_filename,bbox_inches="tight")
    

if __name__=="__main__":
    h5list = glob.glob(output_dir+"*.h5")
    print(h5list)

#    for filename in h5list:
#        print(filename)
#        plot_one(filename)
#        break

#    plot_all(h5list,output_filename=os.path.join(output_dir,"posterior_limits_all.png"))
    plot_all(h5list,output_filename=os.path.join(output_dir,"posterior_limits_q_all.png"))
    
