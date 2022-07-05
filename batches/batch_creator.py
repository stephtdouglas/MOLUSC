# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 11:25:37 2022

@author: Jared
"""
import os
from astropy.table import Table
import numpy as np


csv_path_drive = os.path.expanduser(r'G:\Shared drives\DouglasGroup\Jared Sofair 2022\CSV Files')
csv_path_github = os.path.expanduser(r'C:\Users\Jared\Documents\GitHub\data-parser\CSV Files')
contrast_path = os.path.expanduser(r'G:\Shared drives\DouglasGroup\Jared Sofair 2022\MOLUSC\Data Parser Tables')

anaconda_path = os.path.expanduser(r'C:\Users\Jared\anaconda3')

repo_path = os.path.expanduser(r'C:\Users\Jared\Documents\GitHub\MOLUSC')
gui_path = os.path.join(repo_path, r'code\BinaryStarGUI.py')
batch_path = os.path.join(repo_path, r'batches')

output_path = os.path.expanduser(r'G:\Shared drives\DouglasGroup\Jared Sofair 2022\MOLUSC\MOLUSC Outputs\Tables')

pm = Table.read(os.path.join(csv_path_drive, r'Praesepe_Merged.csv'))
targets = Table.read(os.path.join(csv_path_github, r'targets_abr.csv'))

# Format for running star through MOLUSC in command line
# python [GUI code path] [input data type(s)] [contrast/RV path] [input args] [output path] [other input args]
# Example:
# python "C:/Users/Jared/Documents/GitHub/MOLUSC/code/BinaryStarGUI.py" -v -a --ao "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables/JS355.txt" --filter K --age 0.8000 "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/MOLUSC Outputs/Tables/JS355" 08h40m22.16s -- +18d07m24.8s 3 0.599

def run_batch_stars(stars=["JS355"], analysis_options="ao", write_all=True, extra_output=True, filt=None, age=0.015, companions=3, opsys="win"):
    if opsys=="win":
        # Create .bat file, write line necessary to run the script for Windows
        with open(os.path.join(batch_path, r"batch_runner.bat"), 'w') as f:
            
            # Python seems to have trouble with the back slashes, so I separated this into 3 lines
            f.write(f'call {anaconda_path}')
            f.write(r'\Scripts\activate.bat ')
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
                    f.write(f'--{analysis_options} "{os.path.join(contrast_path, star.replace(" ", "_"))}.txt" --filter {filt} ')
                # if "rv" in analysis_options: # RV
                # if "ruwe" in analysis_options: # RUWE
                # if "gaia" in analysis_options: # Gaia
                
                # Star info (age, output path, ra, dec, # companions, mass)
                # Get ra, dec, and mass
                ra = targets["ra"][np.where(targets["name"] == star)[0][0]]
                dec = targets["de"][np.where(targets["name"] == star)[0][0]]
                mass = np.round(targets["M/Ms"][np.where(targets["name"] == star)[0][0]], 3)
                f.write(f'--age {age} "{os.path.join(output_path, star.replace(" ", "_"))}" {ra} -- +{dec} {companions} {mass}\n\n')
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
        # Create .bat file, write line necessary to run the script for Windows
        with open(os.path.join(batch_path, r"batch_runner.bash"), 'w') as f:
                       
            # Write the command for each star
            for star in stars:
                    
                f.write(f'# {star}\npython "{gui_path}" ') # Label for readability :)
                
                # Write all
                if write_all == True:
                    f.write('-v ')
                # Extra output
                if extra_output == True:
                    f.write('-a ')
                
                # Analysis options and output
                if "ao" in analysis_options: # HRI
                    f.write(f'--{analysis_options} "{os.path.join(contrast_path, star.replace(" ", "_"))}.txt" --filter {filt} ')
                # if "rv" in analysis_options: # RV
                # if "ruwe" in analysis_options: # RUWE
                # if "gaia" in analysis_options: # Gaia
                
                # Star info (age, output path, ra, dec, # companions, mass)
                # Get ra, dec, and mass
                ra = targets["ra"][np.where(targets["name"] == star)[0][0]]
                dec = targets["de"][np.where(targets["name"] == star)[0][0]]
                mass = np.round(targets["M/Ms"][np.where(targets["name"] == star)[0][0]], 3)
                f.write(f'--age {age} "{os.path.join(output_path, star.replace(" ", "_"))}" {ra} -- +{dec} {companions} {mass}\n\n')
        
if __name__ == "__main__": # hehe
    all_targets = targets["name"]
    run_batch_stars(all_targets, age=0.8000, filt="K")