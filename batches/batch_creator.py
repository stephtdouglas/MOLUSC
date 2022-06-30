# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 11:25:37 2022

@author: Jared
"""
import os
from astropy.table import Table


csv_path_drive = os.path.expanduser(r'G:\Shared drives\DouglasGroup\Jared Sofair 2022\CSV Files')
csv_path_github = os.path.expanduser(r'C:\Users\Jared\Documents\GitHub\data-parser\CSV Files')
contrast_path = os.path.expanduser(r'G:\Shared drives\DouglasGroup\Jared Sofair 2022\MOLUSC\Data Parser Tables')

anaconda_path = os.path.expanduser(r'C:\Users\Jared\anaconda3')

repo_path = os.path.expanduser(r'C:\Users\Jared\Documents\GitHub\MOLUSC')
gui_path = os.path.join(repo_path, r'code\BinaryStarGUI.py')
batch_path = os.path.join(repo_path, r'batches')

pm = Table.read(os.path.join(csv_path_drive, r'Praesepe_Merged.csv'))
targets = Table.read(os.path.join(csv_path_github, r'targets_abr.csv'))


# Format for running star through MOLUSC in command line
# python [GUI code path] [input data type(s)] [contrast/RV path] [input args] [output path] [other input args]
# Example:
# python "C:/Users/Jared/Documents/GitHub/MOLUSC/code/BinaryStarGUI.py" -v -a --ao "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables/JS355.txt" --filter K --age 0.8000 "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/MOLUSC Outputs/Tables/JS355" 08h40m22.16s -- +18d07m24.8s 3 0.599

def run_batch_stars(stars=["JS355"], analysis_options="ao", write_all=True, extra_output=True, filt=None, opsys="win"):
    if opsys=="win":
        # Create .bat file, write line necessary to run the script for Windows
        with open(os.path.join(batch_path, r"batch_runner.bat"), 'w') as f:
            f.write(f"call {anaconda_path}\Scripts\activate.bat {anaconda_path}\n\n")
            
            for star in stars:
                
                
                
                # Write the command for each star
                f.write(f"# {star}\n")
                if write_all == True and extra_output == True:
                    f.write(f'python "{gui_path}" -v -a --{analysis_options} "{os.path.join(contrast_path, star)}.txt" --filter {filt} \n\n')
    
    
    
    
    
    
    
    else:
        print(':-(')
        
if __name__ == "__main__":
    all_targets = targets["name"]
    run_batch_stars(all_targets, filt="K")