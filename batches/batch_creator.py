# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 11:25:37 2022

@author: Jared Sofair, Stephanie Douglas
"""
import os, pathlib
from astropy.table import Table
import numpy as np

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
data_path = os.getenv("DATA_PATH")

csv_path_hpc = os.path.join(repo_path, 'csvs')
contrast_path_hpc = os.path.join(repo_path, 'prae_keck')
detection_path_hpc = os.path.join(repo_path, 'prae_keck_detect')
rv_table_path_hpc = os.path.join(repo_path, 'prae_rvs')
slurm_path_hpc = os.path.join(repo_path, 'cluster_scripts')


# repo_path = os.path.expanduser(r'C:/Users/Jared/Documents/GitHub/MOLUSC')
gui_path = os.path.join(repo_path, r'molusc/BinaryStarGUI.py')
batch_path = os.path.join(repo_path, r'batches')
gui_path_hpc = os.path.join(repo_path, r'molusc/BinaryStarGUI.py')
batch_path_hpc = os.path.join(repo_path, r'batches')

output_path_drive = os.path.expanduser('~/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/MOLUSC Outputs/Tables')
# table_output_path_hpc = os.path.join(repo_path, r'../saves/outputs/tables')
# scratch_output_path_hpc = os.path.join(repo_path, r'/scratch/sofairj')
# yml_output_path_hpc = os.path.join(repo_path, r'../saves/outputs/yml')
output_path_hpc = os.path.join(data_path,"molusc_outputs/")

pm = Table.read(os.path.join(csv_path_hpc, r'praesepe_merged.csv'))
targets = Table.read(os.path.join(csv_path_hpc, r'targets_abr.csv'))


# Format for running star through MOLUSC in command line
# python [GUI code path] [input data type(s)] [contrast/RV path] [input args] [output path] [other input args]
# Example:
# python "C:/Users/Jared/Documents/GitHub/MOLUSC/code/BinaryStarGUI.py" -v -a --ao "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/Data Parser Tables/JS355.txt" --filter K --age 0.8000 "G:/Shared drives/DouglasGroup/Jared Sofair 2022/MOLUSC/MOLUSC Outputs/Tables/JS355" 08h40m22.16s -- +18d07m24.8s 3 0.599
ao_path = ''
rv_path = ''


def run_batch_stars(stars=["JS355", "JS364"], yml=True, analysis_options=["ao"], 
                    write_all=True, extra_output=True, filt=None, res=None, 
                    rv_floor=None, companions=10, opsys="linux",
                    comp_file=None):
    """
    Creates a single script to run all stars in serial
    """

    if opsys=="win":
        #%% Create .bat file, write line necessary to run the script for Windows
        with open(os.path.join(batch_path, r"batch_runner.bat"), 'w') as f:
            # pm = Table.read(os.path.join(csv_path_github, r'Praesepe_Merged.csv'))
            # targets = Table.read(os.path.join(csv_path_github, r'targets_abr.csv'))
            
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
                    ao_path = os.path.join(contrast_path, star.replace(" ", "_"))    
                    f.write(f'--ao "{ao_path}.txt" --filter {filt} ')
                if "rv" in analysis_options: # RV
                    rv_path = os.path.join(rv_table_path_hpc, star.replace(" ", "_"))
                    f.write(f'--rv "{rv_path}.txt" --resolution {res} --rv_floor {rv_floor} ')
                if "gaia" in analysis_options: # Gaia
                    f.write('--gaia ')
                
                # Star info (age, output path, ra, dec, # companions, mass)
                # Get ra, dec, and mass
                ra = targets["ra"][np.where(targets["name"] == star)[0][0]]
                dec = targets["de"][np.where(targets["name"] == star)[0][0]]
                mass = np.round(targets["M/Ms"][np.where(targets["name"] == star)[0][0]], 3)
                age = targets["age"][np.where(targets["name"] == star)[0][0]]
                print(age)
                out_path_drive = os.path.join(output_path_drive, star.replace(" ", "_"))
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
            with open(os.path.join(batch_path_hpc, r"batch_runner_yml.bash"), 'w') as f:
                f.write('#!/bin/bash\n\n')
                
                # Write the line for each star
                for star in stars:
                    f.write(f'# {star}\npython {gui_path_hpc} ') # Label for readability :)
                    f.write(f'yml {yml_output_path_hpc}/{star.replace(" ", "_")}_params_output.yml\n\n')
                    
        # You want to create the fml file for the star(s) you intend to run through MOLUSC        
        else:
            with open(os.path.join(slurm_path_hpc, f"batch_runner_cl.bash"), 'w') as f:
                f.write('#!/bin/bash\n\n')
                
                # Write the line for each star
                for star in stars:
                    f.write(f'# {star}\npython {gui_path_hpc} ') # Label for readability :)
                    f.write('cl ')
                
                    # Write all
                    if write_all == True:
                        f.write('-v ')
                    # Extra output
                    if extra_output == True:
                        f.write('-a ')
                    
                    # Analysis options
                    if "ao" in analysis_options: # HRI
                        ao_path = os.path.join(contrast_path_hpc, star.replace(" ", "_"))+".txt"
                        # print(ao_path)
                        if os.path.exists(ao_path):
                            f.write(f'--ao {ao_path} --filter {filt} ')
                    if "rv" in analysis_options: # RV
                        rv_path = os.path.join(rv_table_path_hpc, star.replace(" ", "_"))+".txt"
                        # print(rv_path)
                        if os.path.exists(rv_path):
                            f.write(f'--rv {rv_path} --resolution {res} --rv_floor {rv_floor} ')
                    if "gaia" in analysis_options: # Gaia
                        f.write('--gaia ')

                    if (comp_file is not None) and (os.path.exists(comp_file)):
                        print(comp_file)
                        f.write(f'--comps {comp_file} ')
                    
                    # Star info (age, output path, ra, dec, # companions, mass)
                    # Get ra, dec, and mass
                    ra = targets["ra"][np.where(targets["name"] == star)[0][0]]
                    dec = targets["de"][np.where(targets["name"] == star)[0][0]]
                    mass = np.round(targets["M/Ms"][np.where(targets["name"] == star)[0][0]], 3)
                    age = targets["age"][np.where(targets["name"] == star)[0][0]]
                    
                    out_path_hpc = os.path.join(output_path_hpc, star.replace(" ", "_"))
                    f.write(f'--age {age} -- {out_path_hpc} {ra} +{dec} {companions} {mass}\n\n')


def create_slurm_script(star, yml=True, analysis_options=["ao"], 
                        write_all=True, extra_output=True, filt=None, res=None, 
                        rv_floor=None, companions=10, opsys="linux",
                        comp_file=None):
    """
    Create slurm scripts to run a single star individually
    """
    star_name = star.replace(" ","_")

    # Analysis options
    if "ao" in analysis_options: # HRI
        # If a detection exists, use that instead of the contrast limits
        det_path = os.path.join(detection_path_hpc, star_name)+"_detection.csv"
        # print(ao_path)
        if os.path.exists(det_path):
            ao_path = det_path
        else:
            ao_path = os.path.join(contrast_path_hpc, star_name)+".txt"    
        if os.path.exists(ao_path)==False:
            analysis_options.remove("ao")
    if "rv" in analysis_options: # RV
        rv_path = os.path.join(rv_table_path_hpc, star_name)+".txt"
        # print(rv_path)
        if os.path.exists(rv_path)==False:
            analysis_options.remove("rv")
            # ignoring anything without rvs right now
#            return False


    with open(os.path.join(batch_path_hpc, f"run_{star_name}.sh"), 'w') as f:
        f.write('#!/bin/bash\n#\n')

        f.write(f"#SBATCH --job-name={star_name}\n")
        f.write(f"#SBATCH --output=/data/labs/douglaslab/douglste/script_logs/slurm-%A_{star_name}.out\n")
        f.write("#SBATCH --account=douglste-laf-lab\n")
        f.write("#SBATCH --partition=douglste-laf-lab,compute,unowned\n")
        f.write("#\n")
        f.write("#SBATCH --ntasks=1\n")

        # The configuration depends on what data needs to be analyzed

        cpus = 8
        hours = 0.25

        if companions<=4_000_000:
            cpus = 2
            mem_cpu = 24
            hours += 0.25
        elif companions<=40_000_000:
            cpus = 8
            mem_cpu = 23
            hours += 0.5
        elif companions<=50_000_000:
            print("WARNING: this will require the himem queue, which is currently failing")
            cpus = 6
            mem_cpu = 36
            hours += 0.75
        else:
            print("WARNING: make sure this isn't an overkill request")
            cpus = 6
            mem_cpu = 40
            hours += 3
            
        
        if (comp_file is None) or (os.path.exists(comp_file))==False:
            hours += 0.5

        if "rv" in analysis_options:
            hours += 2.5
            mem_cpu +=30
            cpus +=2
            
        f.write("#SBATCH --cpus-per-task=")
        f.write(f"{cpus}\n")
        f.write("#SBATCH --time=")
        mins = int((hours % 1)*60)
        hours = int(hours)
        f.write(f"{hours}:{mins:02d}:00\n")
        mem_cpu = int(mem_cpu)
        f.write(f"#SBATCH --mem-per-cpu={mem_cpu}gb\n")

        f.write("#SBATCH --mail-type=END,FAIL\n")
        f.write("#SBATCH \n")
        f.write("#SBATCH --mail-user=douglste@lafayette.edu\n\n")

        f.write("srun hostname\n\n")
        f.write("source ~/.bashrc\n\n")
        f.write("source activate molusc\n\n")
#        f.write("cd ~/projects/MOLUSC/\n\n")
        f.write("cd /scratch/ \n\n")

        # make some other things easier by assigning
        # the star name to a bash variable
        f.write(f'sname="{star_name}"\n')
        f.write(f'sdir="douglste_$sname"\n')
        
#        f.write(f"mkdir douglste_{star_name} \n\n")
        f.write('if [ -d "$sdir" ]; then\n')
        f.write('    echo "Deleting old files"\n')
        f.write('    rm -fr "$sdir"\n')
        f.write('fi\n\n')
        f.write('if [ -f "$sname"_RVs.csv ]; then\n')
        f.write('    echo "Deleting old RV file"\n')
        f.write('    rm -f "$sname"_RVs.csv \n')
        f.write('fi\n\n')
        f.write('mkdir "$sdir" \n\n')

        f.write('echo "ls result"\n')
        f.write("ls \n\n")

        comp_file_short = comp_file.split("/")[-1]
        comp_file_scratch = f"/scratch/douglste_{star_name}/{comp_file_short}"
        f.write(f"cp {comp_file} {comp_file_scratch} \n\n")

        
        # Call the python application, set to command line mode
        f.write(f'srun python {gui_path_hpc} ')
        f.write('cl ')
    
        # Write all
        if write_all == True:
            f.write('-v ')
        # Extra output
        if extra_output == True:
            f.write('-a ')
        
        # Analysis options
        if "ao" in analysis_options: # HRI
            # ao_path = os.path.join(contrast_path_hpc, star_name)+".txt"
            # print(ao_path)
            if os.path.exists(ao_path):
                f.write(f'--ao {ao_path} --filter {filt} ')
            if "detect" in ao_path:
                f.write("--detect ")
        if "rv" in analysis_options: # RV
            rv_path = os.path.join(rv_table_path_hpc, star_name)+".txt"
            # print(rv_path)
            if os.path.exists(rv_path):
                f.write(f'--rv {rv_path} --resolution {res} --rv_floor {rv_floor} ')
        if "gaia" in analysis_options: # Gaia
            f.write('--gaia ')

        if (comp_file is not None) and (os.path.exists(comp_file)):
            f.write(f'--comps {comp_file_scratch} ')
        
        # Star info (age, output path, ra, dec, # companions, mass)
        # Get ra, dec, and mass
        ra = targets["ra"][np.where(targets["name"] == star)[0][0]]
        dec = targets["de"][np.where(targets["name"] == star)[0][0]]
        mass = np.round(targets["M/Ms"][np.where(targets["name"] == star)[0][0]], 3)
        age = targets["age"][np.where(targets["name"] == star)[0][0]]
        
#        out_path_hpc = os.path.join(output_path_hpc, star_name)
        out_path_hpc = f"/scratch/douglste_{star_name}/{star_name}"
        f.write(f'--age {age} -- {out_path_hpc} {ra} +{dec} {companions} {mass}\n\n')

        f.write(f"ls /scratch/douglste_{star_name} \n\n")
        f.write(f"mv /scratch/douglste_{star_name}/{star_name}* {output_path_hpc} \n\n")
        f.write(f"rm -f /scratch/douglste_{star_name}/MOLUSC_prior*.hdf5 \n\n")
        f.write(f"rmdir /scratch/douglste_{star_name} -v --ignore-fail-on-non-empty \n\n")

        return True

if __name__ == "__main__": 
    all_targets = targets["name"]
    # rv_targets = 1
    # no_rv_targets = 2
    # run_batch_stars(all_targets, yml=False, analysis_options=["ao"], filt="K", companions=15000, opsys='linux')
#    run_batch_stars(all_targets, yml=False, write_all=True, extra_output=True, 
#                    analysis_options=["ao", "rv", "gaia", "ruwe"], 
#                    filt="K", companions=50_000_000, rv_floor=1000, res=20_000, opsys='linux',
#                    comp_file="/data/douglaslab/douglste/molusc_cache/MOLUSC_prior_v0_Pflat_50M.hdf5")
    # run_batch_stars(stars=["JS355", "JS364"], yml=False, 
    # analysis_options=["ao", "gaia", "ruwe", "rv"], write_all=True, 
    # extra_output=True, filt="K", rv_floor=1000, res=20000, companions=10, opsys="linux")

    with open("submit_all.sh","w") as g:
        for name in all_targets:
            star_name = name.replace(" ","_")
            outfile1 = os.path.join(data_path,f"molusc_outputs/{star_name}_kept.csv")
            outfile2 = os.path.join(data_path,f"molusc_outputs/{star_name}.h5")
            #if os.path.exists(outfile1) or os.path.exists(outfile2):
            #    continue
            
            to_run = create_slurm_script(name, yml=False, write_all=True, extra_output=True, 
                        analysis_options=["ao", "gaia", "ruwe", "rv"], 
                        filt="K", companions=100_000_000, rv_floor=1000, res=20_000, opsys='linux',
                        comp_file=os.path.join(data_path,"molusc_cache/MOLUSC_prior_v0_Pflat_100M.hdf5"))

            if to_run:
                batch_script = os.path.join(batch_path_hpc,f"run_{star_name}.sh")
                g.write(f"sbatch {batch_script}\n")


