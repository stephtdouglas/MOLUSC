import os, sys, glob

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.ascii as at
import corner
import h5py
import yaml

def posterior(post_file,yaml_file,star_name,output_filename):
    
    post = at.read(post_file)
    post["logP(d)"] = np.log10(post["period(days)"])
    post_dict = {}
    for name in ['logP(d)','eccentricity','mass ratio','cos_i']:
        post_dict[name] = np.asarray(post[name])

    fig = corner.corner(post_dict,plot_datapoints=False)

    plt.suptitle(star_name)

    with open(yaml_file, "r") as g:
        data = yaml.safe_load(g)

    ax = fig.axes[2]
    textx = 0.3
    if data["ao_params"]["fit"]==True:
        ax.text(textx,0.8,"AO - Keck")
    if data["rv_params"]["fit"]==True:
        ax.text(textx,0.7,"RV - Hydra")
    if data["gaia_params"]["fit"]==True:
        ax.text(textx,0.6,f"Gaia; RUWE={data['gaia_params']['ruwe']:.1f}")
    ngen = data["num_generated"]
    perc = len(post)/ngen*100
    ax.text(textx,0.5,f"{perc:.1f}% remain of {ngen//1_000_000} Million")
    ax.text(textx,0.4,f"{data['star']['ra']} {data['star']['dec']}")
    ax.text(textx,0.3,f"age: {data['star']['age']:.2f} Gyr")
    ax.text(textx,0.2,f"mass: {data['star']['mass']:.2f} Msun")

    plt.savefig(output_filename,bbox_inches="tight",dpi=600)
    plt.close()

if __name__=="__main__":
    
    output_dir = os.path.expanduser("~/data2/molusc_outputs/")

    files = glob.glob(output_dir+"*kept.csv")

    for file in files:
        star_name = file.split("/")[-1].replace("_kept.csv","")
        output_filename = file.replace(".csv",".png")
        yml_filename = file.replace("_kept.csv","_params_output.yml")
        posterior(file,yml_filename,star_name,output_filename)
