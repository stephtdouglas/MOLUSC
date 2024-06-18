import os

import numpy as np
import yaml

data_dir = os.getenv("MOLOUT", './')

# Convenience functions
def jup_mass_to_sol(jupiter_mass):
    return 9.457e-4*jupiter_mass

def sol_mass_to_jup(solar_mass):
    return solar_mass/9.457e-4

def period_to_a(per):
    # Returns the semi major axis (in AU) for a given period (in days) assuming an equal mass binary with solar mass
    # stars
    G = 39.478  # Gravitational constant in AU^3/years^2*M_solar
    return (np.divide(per, 365)**2 * G * 2/(4*np.pi**2))**(1 / 3)

def a_to_period(a):
    # Returns period (days) of an equal mass, solar-mass binary with a given semi-major axis
    G = 39.478 # Gravitational constant in AU^3/years^2*M_solar
    return np.sqrt((4*np.pi**2 * a**3) / (2*G)) * 365


def get_info(star,output_dir=data_dir):
    output_path1 = output_dir #os.path.join(output_dir, "tables/")
    output_path2 = output_dir #os.path.join(output_dir, "graphs/")
    output_path3 = output_dir #os.path.join(output_dir, "yml/")

    h5file = f'{output_path1}/{star}.h5'
    if os.path.exists(h5file):
        all_file = h5file
        survivors_file = h5file
    else:
        survivors_file = os.path.expanduser(f'{output_path1}/{star}_kept.csv')
        all_file = os.path.expanduser(f'{output_path1}/{star}_all.csv')
    yaml_file = os.path.expanduser(f'{output_path3}/{star}_params_output.yml')

    out_file1 = os.path.expanduser(f'{output_path2}/{star}_output_corner.pdf')
    out_file2 = os.path.expanduser(f'{output_path2}/{star}_output_dtct_lims.pdf')
    out_file3 = os.path.expanduser(f'{output_path2}/{star}_output_srv.pdf')
    
    # Get mass and n  
    with open(yaml_file, 'r') as f:
        try:
            yaml_data = yaml.safe_load(f)
            mass = yaml_data["star"]["mass"]
            n = yaml_data["num_generated"]
        except yaml.YAMLError as exc:
            print(exc)
            
    info = [survivors_file, all_file, out_file1, out_file2, out_file3, mass, n]
    # info[0] = survivors_file
    # info[1] = all_file
    # info[2] = out_file1
    # info[3] = out_file2
    # info[4] = out_file3
    # info[5] = mass
    # info[6] = n
    return info
