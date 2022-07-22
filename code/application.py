from datetime import datetime
import numpy as np
import scipy as scipy
import scipy.stats as stats
from astropy.io import ascii
import argparse
import gc
import warnings
import os
import yaml
from ao import AO
from companions import Companions
from gui import GUI
from ruwe import RUWE
from rv import RV
from astropy.utils.exceptions import AstropyWarning
import logging
logging.basicConfig(format='%(asctime)s %(message)s', encoding='utf-8', level=logging.DEBUG)
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

today = datetime.today().isoformat().split("T")[0]
global repo_path
repo_path = os.getenv('MOLOC').replace("\\", "/")

class Application:
    # Input
    input_args = []
    using_gui = True
    # Class variables
    ao_filename = ['']
    filter = ''
    rv_filename = ''
    resolution = 50000
    prefix = ''
    ruwe_check = False
    gaia_check = False
    # Star parameters
    star_ra = '00h00m00.00s'
    star_dec = '00d00m00.0s'
    star_mass = 0
    star_age = 5
    added_jitter = 20
    rv_floor = 20
    # Code Parameters
    num_generated = 0
    limits = []
    pd_mu = 5.03
    pd_sig = 2.28
    q_exp = 0.0
    gaia_limit = 18.
    extra_output = False
    all_output = False
    # Reject lists
    ao_reject_list = []
    rv_reject_list = []
    jitter_reject_list = []
    ruwe_reject_list = []
    gaia_reject_list = []

    # Class functions
    def __init__(self, input_args):
        self.input_args = input_args

    def start(self):
        # Start has to handle the input arguments
        if len(self.input_args) == 1:
            # run with the gui
            self.using_gui = True
            self.gui = GUI(self)
            self.gui.start()
        else:
            # run from the arguments only
            self.using_gui = False
            if self.parse_input():
                self.run()
    
    def get_inputs(self):
        # Get variables from the GUI
        # Analysis Options
        self.ao_filename = self.gui.get_ao_filename()
        self.rv_filename = self.gui.get_rv_filename()
        self.ruwe_check = self.gui.get_ruwe()
        self.gaia_check = self.gui.get_gaia()
        self.filter = self.gui.get_filter()
        self.resolution = self.gui.get_resolution()
        # Output Options
        self.prefix = self.gui.get_prefix()
        self.extra_output = self.gui.get_extra()
        self.all_output = self.gui.get_all_out_bool()
        # Star Information
        self.star_ra = self.gui.get_ra()
        self.star_dec = self.gui.get_dec()
        self.star_age = self.gui.get_age()
        self.star_mass = self.gui.get_mass()
        self.num_generated = self.gui.get_num_generated()
        self.added_jitter = self.gui.get_added_jitter()
        self.rv_floor = self.gui.get_rv_floor()
        # Limits
        self.limits = self.gui.get_limits()
        # Period distribution
        self.pd_mu, self.pd_sig = self.gui.get_p_dist()
        # Mass Ratio distribution
        self.q_exp = self.gui.get_q_dist()
        # Gaia Completeness
        self.gaia_limit = float(self.gui.get_gaia_limit())
        return
    
    def run(self):
        gc.collect()
        if self.using_gui:
            # Clear display
            self.gui.gui_print('clc')
            # Update status
            self.gui.update_status('Running')
            # Get user inputs from GUI
            self.get_inputs()

        # Print out star info  # TESTING
        self.print_out(('Star RA: ' + self.star_ra))
        self.print_out(('Star DEC: ' + self.star_dec))
        self.print_out(('Star Mass: ' + str(self.star_mass)))

        # Generate Companions
        self.print_out('Generating Companions..')
        comps = Companions(self.num_generated, self.limits, self.star_mass, self.pd_mu, self.pd_sig, self.q_exp)
        failure = self.error_check(comps.generate())
        if failure: return
        self.print_out('Companions Generated')
        self.print_out(('Number of companions: ' + str(self.num_generated)))

        # Decide which parts to run, and run them
        #    AO
        if self.ao_filename[0]:
            ao_reject_lists = []
            for i in range(len(self.ao_filename)):
                if self.ao_filename[i]:
                    self.print_out(f'Analyzing contrast curve in {self.ao_filename[i]}')
                    ao = AO(self.ao_filename[i], comps, self.star_mass, self.star_age, self.star_ra, self.star_dec, self.filter[i])
                    # Determine distance
                    failure = self.error_check(ao.get_distance(self.star_ra, self.star_dec, self.parallax))
                    if failure: return
                    if self.extra_output: self.print_out(('Calculated distance to star: %.0f pc' % (ao.star_distance*4.84e-6)))
                    # Read contrast file
                    failure = self.error_check(ao.read_contrast())
                    if failure: return
                    if self.extra_output:
                        self.print_out('AO Contrast Loaded')
                    # Perform test
                    result = ao.analyze()
                    failure = self.error_check(result)
                    if failure: return
                    ao_reject_lists.append(result)
                    if self.extra_output:
                        self.print_out('Star Model Mag: %.2f' %(ao.star_model_mag))
                        self.print_out(('AO Low Mass Limit: %.3f' %(ao.low_mass_limit)))
            # Combine reject lists from all AO tests into one
            self.ao_reject_list = np.logical_or.reduce(ao_reject_lists)
        else:
            self.ao_reject_list = np.array([False]*self.num_generated)

        #   RV and Jitter
        if self.rv_filename:
            # Run RV (without Jitter)
            self.print_out('Analyzing RV...')
            rv = RV(self.rv_filename, self.resolution, comps, self.star_mass, self.star_age, added_jitter=self.added_jitter, rv_floor=self.rv_floor, extra_output=self.extra_output)
            # Read in the RV file
            failure = self.error_check(rv.read_in_rv())
            if failure: return
            if self.extra_output: self.print_out('RV Measurements Loaded')
            # Run analysis
            self.rv_reject_list = rv.analyze_rv()
            self.jitter_reject_list = np.array([False]*self.num_generated)
        else:
            self.rv_reject_list = np.array([False]*self.num_generated)
            self.jitter_reject_list = np.array([False] * self.num_generated)

        #   RUWE
        if self.ruwe_check:
            self.print_out('Analyzing RUWE...')
            ruwe = RUWE(self.star_ra, self.star_dec, self.star_age, self.star_mass, comps)
            # Read in RUWE distribution and Normalization tables
            failure = self.error_check(ruwe.read_dist())
            if failure: return
            # Get Gaia information
            if np.isfinite(self.gmag):
                ruwe.gmag = self.gmag
                ruwe.color = self.color
                ruwe.n_good_obs = self.n_good_obs
                ruwe.astrometric_chi2 = self.astrometric_chi2
                ruwe.parallax = self.parallax
                ruwe.parallax_error = self.parallax_error
                ruwe.ln_ruwe = self.ln_ruwe
                logging.debug(f"self.ln_ruwe vs. ruwe.ln_ruwe round 1: {self.ln_ruwe} vs. {ruwe.ln_ruwe}")

            failure = self.error_check(ruwe.get_gaia_info())
            if failure: return
            # Perform Test
            self.ruwe_reject_list = ruwe.analyze()
            if self.extra_output:
                self.print_out(('The star has ln(ruwe) of %f.' % (ruwe.ln_ruwe)))
            logging.debug(f"self.ln_ruwe vs. ruwe.ln_ruwe round 2: {self.ln_ruwe} vs. {ruwe.ln_ruwe}")
        else:
            self.ruwe_reject_list = np.array([False]*self.num_generated)
            logging.debug(f"n is not finite!: {self.ln_ruwe} vs. {ruwe.ln_ruwe}")

        #   Gaia Contrast
        if self.gaia_check:
            self.print_out('Analyzing Gaia Contrast...')
            #todo improve gaia contrast
            gaia = AO(f'{os.path.join(repo_path, "code/gaia_contrast.txt")}', comps, self.star_mass, self.star_age, self.star_ra, self.star_dec, 'G', gaia=True)
            # Determine distance
            failure = self.error_check(gaia.get_distance(self.star_ra, self.star_dec, self.parallax))
            if failure: return
            if self.extra_output: self.print_out(('Calculated distance to star: %.0f pc' % (gaia.star_distance*4.84e-6)))
            # Read contrast file
            failure = self.error_check(gaia.read_contrast())
            if failure: return
            if self.extra_output: self.print_out('Gaia Contrast Loaded')
            # Analyze
            self.gaia_reject_list = gaia.analyze_gaia(self.gaia_limit)
        else:
            self.gaia_reject_list = np.array([False]*self.num_generated)

        # Check successes
        w = [True if x == -1 else False for x in self.ao_reject_list]
        x = [True if x == -1 else False for x in self.rv_reject_list]
        y = [True if x == -1 else False for x in self.jitter_reject_list]
        z = [True if x == -1 else False for x in self.ruwe_reject_list]
        if all(w) or all(x) or all(y) or all(z):
            self.gui.update_status('Finished - Unsuccessful')
            if all(w):
                self.gui.gui_print('AO Problem')
            if all(x):
                self.gui.gui_print('RV Problem')
            if all(y):
                self.gui.gui_print('Jitter problem')
            if all(z):
                self.gui.gui_print('RUWE problem')
            return

        # Put together reject lists and output
        keep = np.invert(self.ao_reject_list)*np.invert(self.rv_reject_list)*np.invert(self.jitter_reject_list)*np.invert(self.ruwe_reject_list)*np.invert(self.gaia_reject_list)
        num_kept = sum([1 for x in keep if x])
        self.print_out(('Total number of surviving binaries: ' + str(num_kept)))

        # Write out files
        #   Write out the survivors file
        cols = ['mass ratio', 'period(days)', 'semi-major axis(AU)', 'cos_i', 'eccentricity', 'arg periastron', 'phase']
        keep_table = np.vstack((comps.get_mass_ratio()[keep], comps.get_P()[keep], comps.get_a()[keep],
                                comps.get_cos_i()[keep], comps.get_ecc()[keep], comps.get_arg_peri()[keep],
                                comps.get_phase()[keep]))
        if self.ao_filename[0]:
            cols = cols + ['Projected Separation(AU)','Model Contrast']
            keep_table = np.vstack((keep_table, np.array(ao.pro_sep)[keep], np.array(ao.model_contrast)[keep]))
        elif self.gaia_check:
            cols = cols + ['Projected Separation(AU)', 'DeltaG']
            keep_table = np.vstack((keep_table, np.array(gaia.pro_sep)[keep], np.array(gaia.model_contrast)[keep]))
        if self.ruwe_check:
            if 'Projected Separation(AU)' in cols and 'DeltaG' not in cols:
                cols = cols + ['DeltaG','Predicted RUWE']
                keep_table = np.vstack((keep_table, np.array(ruwe.delta_g)[keep], np.array(ruwe.predicted_ruwe)[keep]))
            elif 'Projected Separation(AU)' in cols and 'DeltaG' in cols:
                cols = cols + ['Predicted RUWE']
                keep_table = np.vstack((keep_table, np.array(ruwe.predicted_ruwe)[keep]))
            else:
                cols = cols + ['Projected Separation(AU)','DeltaG','Predicted RUWE']
                keep_table = np.vstack((keep_table, np.array(ruwe.projected_sep)[keep], np.array(ruwe.delta_g)[keep], np.array(ruwe.predicted_ruwe)[keep]))
        if self.rv_filename:
            # # Write out RV calculations, makes v. large files
            # cols = cols + ['rv'+str(i) for i in range(0, len(rv.MJD))]
            # keep_table = np.vstack((keep_table, np.transpose(np.array(rv.predicted_RV)[keep])))
            cols = cols + ['RV Amplitude','Binary Type']
            keep_table = np.vstack((keep_table, np.array(rv.amp)[keep], np.array(rv.b_type)[keep]))
        keep_table = np.transpose(keep_table)
        
        ascii.write(keep_table, (self.prefix + "_kept.csv"), format='csv', names=cols, overwrite=True)

        self.print_out(('Surviving binary parameters saved to: ' + self.prefix + '_kept.csv'))
        
       

        
        #  Write out the input file
        if self.all_output:
            cols = ['mass ratio', 'period(days)', 'semi-major axis(AU)', 'cos_i', 'eccentricity', 'arg periastron', 'phase']
            all_table = np.vstack((comps.get_mass_ratio(), comps.get_P(), comps.get_a(), comps.get_cos_i(),
                                   comps.get_ecc(), comps.get_arg_peri(), comps.get_phase()))
            if self.ao_filename[0]:
                cols = cols + ['Projected Semparation(AU)', 'Model Contrast'] + [('AO Rejected ' + str(i+1)) for i in range(len(ao_reject_lists))]
                all_table = np.vstack((all_table, ao.pro_sep, ao.model_contrast, ao_reject_lists))
            if self.rv_filename:
                # Write out RV calculations
                cols = cols + ['RV Amplitude','Binary Type','RV Rejected']
                all_table = np.vstack((all_table, rv.amp, rv.b_type, self.rv_reject_list))
                # Uncomment to  include RV Calcualations in output file (makes it obnoxiously large)
                # cols = cols + ['rv' + str(i) for i in range(0, len(rv.MJD))]
                # all_table = np.vstack((all_table, np.transpose(np.array(rv.predicted_RV))))
            if self.ruwe_check:
                if 'Projected Separation(AU)' not in cols:
                    cols = cols + ['Projected Separation(AU)','DeltaG','Predicted RUWE','RUWE Rejected','RUWE Rejection Prob']
                    all_table = np.vstack((all_table, ruwe.projected_sep, ruwe.delta_g, ruwe.predicted_ruwe, self.ruwe_reject_list, ruwe.rejection_prob))
                else:
                    cols = cols + ['DeltaG','Predicted RUWE','RUWE Rejected','RUWE Rejection Prob']
                    all_table = np.vstack((all_table, ruwe.delta_g, ruwe.predicted_ruwe, self.ruwe_reject_list, ruwe.rejection_prob))

            if self.gaia_check:
                if 'Model Contrast' in cols:
                    # The AO test has been run, projected sep and contrast already included
                    cols = cols + ['Gaia Rejected']
                    all_table = np.vstack((all_table, self.gaia_reject_list))
                elif 'DeltaG' in cols or 'Projected Separation(AU)'in cols:
                    # RUWE test has been run, projected sep and contrast already included
                    cols = cols + ['Gaia Rejected']
                    all_table = np.vstack((all_table, self.gaia_reject_list))
                else:
                    cols = cols + ['Projected Separation(AU)', 'DeltaG', 'Gaia Rejected']
                    all_table = np.vstack((all_table, gaia.pro_sep, gaia.model_contrast, self.gaia_reject_list))
                    
            # Add column containing full rejection info
            cols = cols + ['Full Rejected']
            all_table = np.vstack((all_table, np.invert(keep)))
            all_table = np.transpose(all_table)
                
            ascii.write(all_table, (self.prefix + "_all.csv"), format='csv', names=cols, overwrite=True)


            self.print_out(('Generated binary parameters saved to: ' + self.prefix + '_all.csv'))

        # Write out run inputs to a yaml file
        if self.limits[3] == "transit":
            is_transit = True
        else:
            is_transit = False
            
        # Check if you are using the AO datatype
        if len(self.ao_filename[0]) > 3:
            is_ao = True
        else:
            is_ao = False
        
        # Check if you are using the RV datatype
        if (type(self.rv_filename) == bool) and (self.rv_filename == True):
            is_rv = True
        elif (type(self.rv_filename) == bool) and (self.rv_filename == False):
            is_rv = False
        elif (type(self.rv_filename) != bool) and (len(self.rv_filename) > 2):
            is_rv = True
        else:
            is_rv = False


        yaml_data = {"run_date": today,
                     "file_prefix": self.prefix,
                     "num_generated": self.num_generated,
                     "transit": is_transit,

                     "star":{"ra": self.star_ra,
                             "dec": self.star_dec,
                             "mass":self.star_mass,
                             "age": self.star_age,
                             "jitter": self.added_jitter},
                     
                     "rv_params":{"fit": is_rv,
                                  "rv_file": self.rv_filename,
                                  "rv_floor": self.rv_floor,
                                  "resolution": self.resolution},
                     
                     "ao_params":{"fit": is_ao,
                                  "ao_file": self.ao_filename,
                                  "filter": self.filter},
                     
                     "ruwe_params":{"fit": self.ruwe_check},

                     
                     "gaia_params":{"fit": self.gaia_check,
                                    "id": self.gaia_id,
                                    "gmag": self.gmag,
                                     "color": self.color,
                                    "n_good_obs": self.n_good_obs,
                                    "astrometric_chi2": self.astrometric_chi2,
                                    "parallax": self.parallax,
                                    "parallax_error": self.parallax_error,
                                    "ruwe": float(np.exp(self.ln_ruwe)),
                                    "ln_ruwe": float(self.ln_ruwe),
                                    }
                    }
        logging.debug(f"\nyaml_data: {yaml_data}\n")
        
        yaml_path = os.path.expanduser(self.prefix).replace("tables", "yml").replace(r"\\", "/")
        yaml_fname = f"{yaml_path}_params_output.yml"
        with open(yaml_fname,"w") as f:
            yaml.dump(yaml_data,f)
        self.print_out(('Run parameters saved to: ' + self.prefix + '_params_output.yml'))


        if self.using_gui: self.gui.update_status('Finished - Successful')
        self.restore_defaults()
        return


    def error_check(self, error_code):
        # If the return value is a list, array, zero or none, no error was found and the process completed successfully
        if isinstance(error_code, list) or isinstance(error_code, np.ndarray):
            return False
        if error_code == 0 or error_code is None:
            return False
        # todo: Add more clear commenting
        # todo: Add error message for cases of too young or too old in AO
        elif error_code == -1:
            self.print_out('ERROR: File Not Found\nCheck filename and try again')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsucessful')
            return True
        elif error_code == -2:
            self.print_out('ERROR: RV file in incorrect format')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsuccessful')
            return True
        elif error_code == -11:
            self.print_out('Period, mass ratio and semi-major axis cannot all be fixed at the same time.')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsuccessful')
            return True
        elif error_code == -12:
            self.print_out('Incompatible limits on period and semi-major axis.')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsuccessful')
            return True
        elif error_code == -21:
            # Contrast file not found
            self.print_out('ERROR: Contrast file Not Found\nCheck filename and try again')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsucessful')
            return True
        elif error_code == -22:
            self.print_out('ERROR: AO file in incorrect format')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsuccessful')
            return True
        elif error_code == -23:
            # AO, mass higher than model grid
            self.print_out('ERROR: The mass of the primary star is larger than the model grid for that age can handle')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsuccessful')
            return True
        elif error_code == -24:
            # AO, mass lower than model grid
            self.print_out('ERROR: The mass of the primary star is smaller than the model grid for that age can handle')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsuccessful')
            return True
        elif error_code == -31:
            # RV file not found
            self.print_out('ERROR: RV file Not Found\nCheck filename and try again')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsucessful')
            return True
        elif error_code == -32:
            self.print_out('ERROR: RV file in incorrect format')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsuccessful')
            return True
        elif error_code == -41:
            # RUWE Distribution file not found
            self.print_out("""ERROR: RUWE Normalization file not found\nThe file should be named table_u0_g_col.txt 
            and located in the folder the code is being run from.""")
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsucessful')
            return True
        elif error_code == -42:
            # RUWE Norm file wrong format
            self.print_out("""ERROR: RUWE Normalization file is wrong format.""")
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsucessful')
            return True
        elif error_code == -43:
            # RUWE Distribution file not found
            self.print_out("""ERROR: RUWE distribution file not found\nThe file should be named RuweTableGP.txt 
            and located in the folder the code is being run from.""")
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsucessful')
            return True
        elif error_code == -44:
            # RUWE Dist file wrong format
            self.print_out("""ERROR: RUWE distribution file is wrong format.""")
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsucessful')
            return True
        elif error_code == -51:
            self.print_out('Unable to find source matching coordinates in Gaia')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsuccessful')
            return True
        elif error_code == -52:
            # This error code indicates a warning should be printed, but does not necessitate failure
            self.print_out('WARNING: More than one source within 10 arcseconds of given coordinates. Closest one used.')
            return False
        elif error_code == -53:
            self.print_out('Unable to calculate RUWE for given star. The magnitude or color is outside of bounds')
            if self.using_gui:
                self.gui.update_status('Finished - Unsuccessful')
            else:
                self.print_out('Finished - Unsuccessful')
            return True
        return

    def parse_input(self):
        # Read in command line inputs
        # Create Parser
        parser = argparse.ArgumentParser(description='Find limits on the orbital parameters of unseen companions')

        subparsers = parser.add_subparsers(dest='command')
        my_parser = subparsers.add_parser('cl')
        yaml_parser = subparsers.add_parser('yml')

        #### One option: inputting all info on the command line
        # Add arguments
        #  Required arguments (Positional)
        my_parser.add_argument('prefix', help='The prefix to be used on the output files')
        my_parser.add_argument('ra', help='The RA of the star in hms format', metavar='ra[hms]')
        my_parser.add_argument('dec', help='The Dec of the star in dms format', metavar='dec[dms]')
        my_parser.add_argument('n', help='The number of companions to generate', type=int)
        my_parser.add_argument('mass', help='The mass of the primary in solar masses', type=float, metavar='mass[M_sun]')
        # Optional arguments
        my_parser.add_argument('--age', help='The age of the star in Gyr', type=float, required=False, default=5,
                               metavar='AGE[Gyr]')
        my_parser.add_argument('--jitter', help='The RV jitter to be added in quadrature to the error in m/s',
                               required=False, default=20, metavar='JITTER[m/s]')
        my_parser.add_argument('--rv_floor', help='The lowest RV semi-amplitude which can be rejected in m/s',
                               required=False, default=20, metavar='RV FLOOR[m/s]')
        # Analysis Options
        my_parser.add_argument('--rv', help='The path to the file containing the RV data', required=False,
                               metavar='RV_PATH', default='')
        my_parser.add_argument('--resolution', help='The spectral resolution of the RV data', required=False,
                               metavar='RV_PATH', default=50000)
        my_parser.add_argument('--ao', nargs="*", type=str, help='The path to the file(s) containing the AO data', required=False,
                               metavar='AO_PATH', default=[''])
        my_parser.add_argument('--filter', nargs="*", type=str, help='The filter in which the AO data was taken', required=False,
                               choices=['J','K','H','G', 'Bp','Rp','R','I','L','LL','M'])
        my_parser.add_argument('--ruwe', action='store_true', help='Apply the RUWE test')
        my_parser.add_argument('--gaia', action='store_true', help='Apply the GAIA contrast test')
        # Other options
        my_parser.add_argument('-v','--verbose', action='store_true', help='Turn on extra output')
        my_parser.add_argument('-a', '--all', action='store_true', help='Write out all generated companions')
        my_parser.add_argument('--transit', action='store_true', help='Turn on transit limits')
        
        #### Second option: provide a yaml file with all options
        yaml_parser.add_argument("yml_file",help="full path to the yaml file with input values",
                           type=str)
        yaml_parser.add_argument('-v','--verbose', action='store_true', help='Turn on extra output')
        yaml_parser.add_argument('-a', '--all', action='store_true', help='Write out all generated companions')
        
        # Run parser
        args = parser.parse_args()
        
        if args.command == "cl":
            # Input the inputs
            #  Analysis Options
            self.rv_filename = args.rv
            self.resolution = float(args.resolution)
            self.ao_filename = args.ao
            # logging.info(f"ao_filename full list: {self.ao_filename}")
            # logging.info(f"First item in ao_filename list: {self.ao_filename[0]}")
            # logging.info(f"Length of ao_filename: {len(self.ao_filename)}")
            self.filter = args.filter
            # logging.info(f"args.filter: {args.filter}")
            # logging.info(f"First item in args.filter: {args.filter[0]}")
            self.ruwe_check = args.ruwe
            self.gaia_check = args.gaia
            # Output Options
            self.prefix = args.prefix
            self.extra_output = args.verbose
            self.all_output = args.all
            # Stellar Info
            self.num_generated = args.n
            self.star_ra = args.ra
            self.star_dec = args.dec
            self.star_mass = args.mass
            self.star_age = args.age
            self.added_jitter = float(args.jitter)
            self.rv_floor = args.rv_floor
            # Limits
            self.limits = [None]*21
            if args.transit:
                self.limits[3] = 'transit'
            
            # Gaia data
            self.gaia_id = -99
            self.gmag = np.nan
            self.color = np.nan
            self.n_good_obs = 0
            self.astrometric_chi2 = np.nan
            self.parallax = np.nan
            self.parallax_error = np.nan
            self.ln_ruwe = np.nan
            logging.info(f"\n----------------------------------------------ln_ruwe v2: {self.ln_ruwe}\n")

        elif args.command == "yml":

            with open(args.yml_file,"r") as f:
                data = yaml.safe_load(f)
#             print(data)
#             print(f"This is the jitter::::::::::::::::::::: {data['Star']}")
            # Input the inputs
            #  Analysis Options
            if data["rv_params"]["fit"]==True:
                self.rv_filename = data["rv_params"]["rv_file"]
                self.resolution = data["rv_params"]["resolution"]
                self.rv_floor = data["rv_params"]["rv_floor"]
            else:
                self.rv_filename = False
                self.resolution = 50000
                self.rv_floor = 20

            if data["ao_params"]["fit"]==True:
                #if len(self.ao_filename) > 1:
                self.ao_filename = data["ao_params"]["ao_file"]
                #else:
                    #self.ao_filename = data["ao_params"]["ao_file"]

                self.filter = data["ao_params"]["filter"]
            else:
                self.ao_filename = [False]
                self.filter = None

            self.ruwe_check = data["ruwe_params"]["fit"]
            self.gaia_check = data["gaia_params"]["fit"]

            # Gaia data
            if np.isfinite(data["gaia_params"]["gmag"]):
                self.gmag = data["gaia_params"]["gmag"]
                self.color = data["gaia_params"]["color"]
                self.n_good_obs = data["gaia_params"]["n_good_obs"]
                self.astrometric_chi2 = data["gaia_params"]["astrometric_chi2"]
                self.parallax = data["gaia_params"]["parallax"]
                self.parallax_error = data["gaia_params"]["parallax_error"]
                self.ln_ruwe = np.log(data["gaia_params"]["ruwe"])
                self.gaia_id = data["gaia_params"]["id"]
                
                #logging.debug(f"\n----------------------------------------------ln_ruwe: {self.ln_ruwe}\n")
            else:
                self.gaia_id = -99
                self.gmag = np.nan
                #logging.debug(f"\n----------------------------------------------ln_ruwe: {self.ln_ruwe}\n")

                self.color = np.nan
                self.n_good_obs = 0
                self.astrometric_chi2 = np.nan
                self.parallax = np.nan
                self.parallax_error = np.nan
                self.ln_ruwe = np.nan

            # Output Options
            self.prefix = data["file_prefix"]
            self.extra_output = args.verbose
            self.all_output = args.all
            # Stellar Info
            self.num_generated = data["num_generated"]
            self.star_ra = data["star"]["ra"]
            self.star_dec = data["star"]["dec"]
            self.star_mass = data["star"]["mass"]
            self.star_age = data["star"]["age"]
            self.added_jitter = data["star"]["jitter"]
            # Limits
            self.limits = [None]*21
            if data["transit"]==True:
                self.limits[3] = 'transit'
        else:
            print('Invalid command')
            return False

        print(self.ruwe_check,self.gaia_check)


        if not self.ao_filename[0] and not self.rv_filename and not self.ruwe_check and not self.gaia_check and not self.added_jitter:
            print('At least one analysis type needs to be chosen')
            return False
        if self.ao_filename[0] and not self.filter:
            print('AO given without filter')
            return False
        if self.ao_filename[0] and not (self.filter[0] in ['J','K','H','G', 'Bp','Rp','R','I','L','LL','M']):
            print(self.filter[0])
            print(self.filter[0] in ['J','K','H','G', 'Bp','Rp','R','I','L','LL','M'])
            print("Unknown filter. Filter must be one of ['J','K','H','G', 'Bp','Rp','R','I','L','LL','M']")
            return False

        return True

    def print_out(self, message):
        if self.using_gui:
            self.gui.gui_print(message)
        else:
            print(message)
        return

    def restore_defaults(self):
        # Input
        self.input_args = []
        self.using_gui = True
        # Class variables
        self.ao_filename = ''
        self.filter = ''
        self.rv_filename = ''
        self.prefix = ''
        self.ruwe_check = False
        self.gaia_check = False
        # Star parameters
        self.star_ra = '00h00m00.00s'
        self.star_dec = '00d00m00.0s'
        self.star_mass = 0
        self.star_age = 5
        self.added_jitter = 20
        self.rv_floor = 20
        # Code Parameters
        self.num_generated = 0
        self.limits = []
        self.extra_output = False
        self.all_output = False
        # Reject lists
        self.ao_reject_list = []
        self.rv_reject_list = []
        self.jitter_reject_list = []
        self.ruwe_reject_list = []
        self.gaia_reject_list = []
        return