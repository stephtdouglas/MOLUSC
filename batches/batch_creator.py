# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 11:25:37 2022

@author: Jared
"""
import os
from astropy.table import Table
import numpy as np


csv_path_drive = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/CSV Files')
csv_path_github = os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/data-parser/CSV Files')
contrast_path = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables')
csv_path_hpc = os.path.expanduser(r'csvs')
contrast_path_hpc = os.path.expanduser(r'contrasts')

anaconda_path = os.path.expanduser(r'C:/Users/Jared/anaconda3')

repo_path = os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/MOLUSC')
gui_path = os.path.join(repo_path, r'code/BinaryStarGUI.py').replace("\\", "/")
batch_path = os.path.join(repo_path, r'batches').replace("\\", "/")
repo_path_hpc = os.path.expanduser(r'.')
gui_path_hpc = os.path.join(repo_path_hpc, r'code/BinaryStarGUI.py').replace("\\", "/")
batch_path_hpc = os.path.join(repo_path_hpc, r'batches').replace("\\", "/")

output_path_drive = os.path.expanduser(r'G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/MOLUSC Outputs/Tables')
output_path_hpc = os.path.expanduser(r'../saves/outputs/tables')

pm = Table.read(os.path.join(csv_path_hpc, r'praesepe_merged.csv').replace("\\", "/"))
targets = Table.read(os.path.join(csv_path_hpc, r'targets_abr.csv').replace("\\", "/"))


# Format for running star through MOLUSC in command line
# python [GUI code path] [input data type(s)] [contrast/RV path] [input args] [output path] [other input args]
# Example:
# python "C:/Users/Jared/Documents/GitHub/MOLUSC/code/BinaryStarGUI.py" -v -a --ao "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables/JS355.txt" --filter K --age 0.8000 "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/MOLUSC Outputs/Tables/JS355" 08h40m22.16s -- +18d07m24.8s 3 0.599

def run_batch_stars(stars=["JS355"], yml=True, analysis_options=["ao"], write_all=True, extra_output=True, filt=None, companions=3, opsys="win"):
    if opsys=="win":
        #%% Create .bat file, write line necessary to run the script for Windows
        with open(os.path.join(batch_path, r"batch_runner.bat").replace("\\", "/"), 'w') as f:
            # pm = Table.read(os.path.join(csv_path_github, r'Praesepe_Merged.csv').replace("\\", "/"))
            # targets = Table.read(os.path.join(csv_path_github, r'targets_abr.csv').replace("\\", "/"))
            # Python seems to have trouble with the back slashes, so I separated this into 3 lines
            f.write(f'call {anaconda_path}')
            f.write(r'/Scripts/activate.bat ')
            f.write(f'{anaconda_path}\n\n')
            
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
                # if "rv" in analysis_options: # RV
                #   f.write(r'--rv "{os.path.join(rv_path, star.replace(" ", "_")).replace("\\", "/")}.txt" --resolution 50000')
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
        #%% Create .bat file, write line necessary to run the script for Linux
        with open(os.path.join(batch_path_hpc, r"batch_runner.bash").replace("\\", "/"), 'w') as f:
            f.write('#!/bin/bash\n')
            # Write the command for each star
            for star in stars:
                f.write(f'# {star}\npython "{gui_path_hpc}" ') # Label for readability :)
                
                if yml:
                    f.write(f'yml "../saves/outputs/yml/{star}_params_output.yml"')
                    
                    
                else:
                    f.write('cl ')
                
                    # Write all
                    if write_all == True:
                        f.write('-v ')
                    # Extra output
                    if extra_output == True:
                        f.write('-a ')
                    
                    if "ao" in analysis_options: # HRI
                        ao_path = os.path.join(contrast_path_hpc, star.replace(" ", "_")).replace("\\", "/")
                        f.write(f'--ao "{ao_path}.txt" --filter {filt} ')
                    # if "rv" in analysis_options: # RV
                    #   f.write(r'--rv "{os.path.join(rv_path, star.replace(" ", "_")).replace(" ", "_"))}.txt" --resolution 50000')
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
                    
                    out_path_hpc = os.path.join(output_path_hpc, star.replace(" ", "_")).replace("\\", "/")
                    f.write(f'--age {age} "{out_path_hpc}" {ra} -- +{dec} {companions} {mass}\n\n')
                    
if __name__ == "__main__": # hehe
    all_targets = targets["name"]
    run_batch_stars(all_targets, yml=False, analysis_options=["ao"], filt="K", companions=15000, opsys='linux')