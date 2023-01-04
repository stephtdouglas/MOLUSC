# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 11:25:37 2022

@author: Jared
"""
import os
from astropy.table import Table
import numpy as np

repo_path = os.path.join(os.getenv('MOLOC')).replace("\\", "/")

csv_path_drive = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/CSV Files')
csv_path_github = os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files')
contrast_path = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables')
csv_path_hpc = os.path.join(repo_path, r'csvs').replace("\\", "/")
contrast_path_hpc = os.path.join(repo_path, r'contrasts').replace("\\", "/")
rv_table_path_hpc = os.path.join(repo_path, r'rvs').replace("\\", "/")

anaconda_path = os.path.expanduser(r'C:/Users/Jared/anaconda3')



# repo_path = os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/MOLUSC')
gui_path = os.path.join(repo_path, r'code/BinaryStarGUI.py').replace("\\", "/")
batch_path = os.path.join(repo_path, r'batches').replace("\\", "/")
gui_path_hpc = os.path.join(repo_path, r'code/BinaryStarGUI.py').replace("\\", "/")
batch_path_hpc = os.path.join(repo_path, r'batches').replace("\\", "/")

output_path_drive = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/MOLUSC Outputs/Tables')
table_output_path_hpc = os.path.join(repo_path, r'../saves/outputs/tables').replace("\\", "/")
scratch_output_path_hpc = os.path.join(repo_path, r'/scratch/sofairj')
yml_output_path_hpc = os.path.join(repo_path, r'../saves/outputs/yml').replace("\\", "/")

pm = Table.read(os.path.join(csv_path_hpc, r'praesepe_merged.csv').replace("\\", "/"))
targets = Table.read(os.path.join(csv_path_hpc, r'targets_abr.csv').replace("\\", "/"))


# Format for running star through MOLUSC in command line
# python [GUI code path] [input data type(s)] [contrast/RV path] [input args] [output path] [other input args]
# Example:
# python "C:/Users/Jared/Documents/GitHub/MOLUSC/code/BinaryStarGUI.py" -v -a --ao "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables/JS355.txt" --filter K --age 0.8000 "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/MOLUSC Outputs/Tables/JS355" 08h40m22.16s -- +18d07m24.8s 3 0.599
ao_path = ''
rv_path = ''


def run_batch_stars(stars=["JS355", "JS364"], yml=True, analysis_options=["ao"], write_all=True, extra_output=True, filt=None, res=None, rv_floor=None, companions=10, opsys="linux"):
    if opsys=="win":
        #%% Create .bat file, write line necessary to run the script for Windows
        with open(os.path.join(batch_path, r"batch_runner.bat").replace("\\", "/"), 'w') as f:
            # pm = Table.read(os.path.join(csv_path_github, r'Praesepe_Merged.csv').replace("\\", "/"))
            # targets = Table.read(os.path.join(csv_path_github, r'targets_abr.csv').replace("\\", "/"))
            
            f.write(f'call {anaconda_path} /Scripts/activate.bat {anaconda_path}\n\n')
            
            # Write the command for each star
            for star in stars:
                    
                f.write(f':: {star}\npython "{gui_path}" ') # Label for readability :)
                
                # Write all
                if write_all == True:
                    f.write('-v ')
                # Extra output
                if extra_output == True:
                    f.write('-a ')
                
                # Analysis options and output
                if "ao" in analysis_options: # HRI
                    ao_path = os.path.join(contrast_path, star.replace(" ", "_")).replace("\\", "/")    
                    f.write(f'--ao "{ao_path}.txt" --filter {filt} ')
                if "rv" in analysis_options: # RV
                    rv_path = os.path.join(rv_table_path_hpc, star.replace(" ", "_")).replace("\\", "/")
                    f.write(f'--rv "{rv_path}.txt" --resolution {res} --rv_floor {rv_floor} ')
                if "ruwe" in analysis_options: # RUWE
                    f.write('--gaia ')
                if "gaia" in analysis_options: # Gaia
                    f.write('--ruwe ')
                
                # Star info (age, output path, ra, dec, # companions, mass)
                # Get ra, dec, and mass
                ra = targets["ra"][np.where(targets["name"] == star)[0][0]]
                dec = targets["de"][np.where(targets["name"] == star)[0][0]]
                mass = np.round(targets["M/Ms"][np.where(targets["name"] == star)[0][0]], 3)
                age = targets["age"][np.where(targets["name"] == star)[0][0]]
                print(age)
                out_path_drive = os.path.join(output_path_drive, star.replace(" ", "_")).replace("\\", "/")
                f.write(f'--age {age} "{out_path_drive}" {ra} -- +{dec} {companions} {mass}\n\n')
            f.write('pause')
    # Choices:
        # Analysis options: HRI, RV, RUWE, Gaia
            # HRI:
                # Filter
            # RV:
                # RV file
                # Resolution
                # Transit
            # RUWE:
                # Check RUWE
            # Gaia:
                # Check Gaia
        # Star info
            # Generated companions*
            # RA*
            # DEC*
            # Mass*
            # Age
            # Added jitter
            # RV floor
        # Output
            # File name
            # Write all
            # Extra output
    
    else:
        #%% Create .bat file
        
        # Write the command for each star. This creates different files for creating the yml files and running stars
        # using the yml files
        
        # You have the yml file for the star(s) you intend to run through MOLUSC
        if yml:
            with open(os.path.join(batch_path_hpc, r"batch_runner_yml.bash").replace("\\", "/"), 'w') as f:
                f.write('#!/bin/bash\n\n')
                
                # Write the line for each star
                for star in stars:
                    f.write(f'# {star}\npython "{gui_path_hpc}" ') # Label for readability :)
                    f.write(f'yml "{yml_output_path_hpc}/{star.replace(" ", "_")}_params_output.yml"\n\n')
                    
        # You want to create the fml file for the star(s) you intend to run through MOLUSC        
        else:
            with open(os.path.join(batch_path_hpc, r"batch_runner_cl.bash").replace("\\", "/"), 'w') as f:
                f.write('#!/bin/bash\n\n')
                
                # Write the line for each star
                for star in stars:
                    f.write(f'# {star}\npython "{gui_path_hpc}" ') # Label for readability :)
                    f.write('cl ')
                
                    # Write all
                    if write_all == True:
                        f.write('-v ')
                    # Extra output
                    if extra_output == True:
                        f.write('-a ')
                    
                    # Analysis options
                    if "ao" in analysis_options: # HRI
                        ao_path = os.path.join(contrast_path_hpc, star.replace(" ", "_")).replace("\\", "/")
                        f.write(f'--ao "{ao_path}.txt" --filter {filt} ')
                    if "rv" in analysis_options: # RV
                        rv_path = os.path.join(rv_table_path_hpc, star.replace(" ", "_")).replace("\\", "/")
                        f.write(f'--rv "{rv_path}.txt" --resolution {res} --rv_floor {rv_floor} ')
                    if "gaia" in analysis_options: # Gaia
                        f.write('--gaia ')
                    if "RUWE" in analysis_options: # RUWE
                        f.write('--ruwe ')
                    
                    # Star info (age, output path, ra, dec, # companions, mass)
                    # Get ra, dec, and mass
                    ra = targets["ra"][np.where(targets["name"] == star)[0][0]]
                    dec = targets["de"][np.where(targets["name"] == star)[0][0]]
                    mass = np.round(targets["M/Ms"][np.where(targets["name"] == star)[0][0]], 3)
                    age = targets["age"][np.where(targets["name"] == star)[0][0]]
                    
                    out_path_hpc = os.path.join(scratch_output_path_hpc, star.replace(" ", "_")).replace("\\", "/")
                    f.write(f'--age {age} "{out_path_hpc}" {ra} -- +{dec} {companions} {mass}\n\n')
                        
if __name__ == "__main__": # hehe
    all_targets = targets["name"]
    # rv_targets = 1
    # no_rv_targets = 2
    # run_batch_stars(all_targets, yml=False, analysis_options=["ao"], filt="K", companions=15000, opsys='linux')
    # run_batch_stars(all_targets, yml=False, write_all=True, extra_output=True, analysis_options=["ao", "rv", "gaia", "ruwe"], filt="K", companions=10, res=1000, opsys='linux')
    run_batch_stars(stars=["JS355", "JS364"], yml=False, analysis_options=["ao", "gaia", "ruwe", "rv"], write_all=True, extra_output=True, filt="K", rv_floor=1000, res=20000, companions=10, opsys="linux")


