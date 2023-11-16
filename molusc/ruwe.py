from datetime import datetime as dt
import datetime
import warnings
import os, pathlib
import logging

import numpy as np
import scipy as scipy
import scipy.stats as stats
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
from astroquery.gaia import Gaia
from astropy.utils.exceptions import AstropyWarning
from multiprocessing import Process
from multiprocessing.pool import Pool
import multiprocessing as mp
from ao import get_pro_sep
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

today = dt.today().isoformat().split("T")[0]

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent


class RUWE:
    # class variables
    reject_list = []
    star_ra = ''
    star_dec = ''
    star_age = 5
    star_mass = 1
    num_generated = 0
    u0_dictionary = {}
    ln_ruwe = 0.
    binary_prob = 0.
    projected_sep = []

    # class functions
    def __init__(self, star_ra, star_dec, age, mass, companions):
        self.star_ra = star_ra
        self.star_dec = star_dec
        self.star_age = age
        self.star_mass = mass
        self.num_generated = companions.num_generated # originally companions.get_num()
        self.a = companions.a
        self.period = companions.a
        self.e = companions.ecc
        self.phase = companions.phase
        self.arg_peri = companions.arg_peri
        self.cos_i = companions.cos_i
        self.mass_ratio = companions.mass_ratio
        self.gmag = np.nan

    def analyze(self):
        # a. Calculate projected separation
        pro_sep = np.zeros(self.num_generated)
        T_0 = 2457388.5  # epoch 2016.0 in JD, corresponding to Gaia eDR3
        T_0_array = np.full(self.num_generated, fill_value=T_0)
        
        # Parallelization
        # Get projected separation
        star_params = np.array([T_0_array, self.period, self.phase, self.e, self.arg_peri, self.cos_i, self.a])
        all_stars = star_params.T
        
        try:
            cpu_ct = len(os.sched_getaffinity(0))-1
            print(f"Current time: {datetime.datetime.now()} -- RUWE cpu_count HPC:", cpu_ct)
        except AttributeError:
            cpu_ct = mp.cpu_count()-1
            print(f"Current time: {datetime.datetime.now()} -- RUWE cpu_count PC:", cpu_ct)
            
        divisor = int(np.ceil(min(self.num_generated / cpu_ct, 200000)))
        
        with Pool(cpu_ct) as pool:
            pro_sep = pool.starmap(get_pro_sep, all_stars, chunksize=divisor)
            
        self.projected_sep = pro_sep
        # print(self.projected_sep)
        # End parallelization

        # b. Calculate distance
        star_distance = 1 / (self.parallax)  # kpc

        # c. Get Barraffe models, and set up the  interpolation to get magnitude
        model = self.load_stellar_model(self.star_age)
        f_mag =  scipy.interpolate.interp1d(model['M/Ms'], model['G'], fill_value='extrapolate')

        # d. Calculate contrast
        primary_mag = f_mag(self.star_mass)
        companion_mag = np.array([f_mag(self.star_mass*x) for x in self.mass_ratio])
        delta_g = np.subtract(companion_mag, primary_mag)
        self.delta_g = delta_g

        # e. Get ruwe
        ruwe = np.exp(self.ln_ruwe)
        log_ruwe = np.log10(ruwe)

        # f. Get predicted RUWE
        #    convert from mas to AU
        star_distance = star_distance * 2.063e+8 # AU
        self.ruwe_dist['Sep(AU)'] = [star_distance * np.tan(np.radians(x/(3.6e6))) for x in np.power(10, self.ruwe_dist['log(sep)']) ]

        #  2D interpolation functions for ruwe and sigma_ruwe
        x_edges = np.unique(np.array(self.ruwe_dist['Sep(AU)']))
        y_edges = np.unique(np.array(self.ruwe_dist['DeltaG']))
        z = np.reshape(np.array(self.ruwe_dist['log(RUWE)']), [len(y_edges), len(x_edges)])
        z_sigma = np.reshape(np.array(self.ruwe_dist['sigma_log(RUWE)']), [len(y_edges), len(x_edges)])

        f_ruwe = scipy.interpolate.interp2d(x_edges, y_edges, z)
        f_sigma = scipy.interpolate.interp2d(x_edges, y_edges, z_sigma)
        
        logging.info(f"projected seps round 2: {self.projected_sep}")
        logging.info(f"delta g: {self.delta_g}")
        logging.info(f"test: {[f_ruwe(self.projected_sep[i], delta_g[i]) for i in range(self.num_generated)]}")

        # pred_log_ruwe = f_ruwe(self.projected_sep, delta_g)
        # self.predicted_ruwe = 10**pred_log_ruwe
        # pred_sigma = f_sigma(self.projected_sep, delta_g)
        
        pred_log_ruwe = np.concatenate([f_ruwe(self.projected_sep[i], delta_g[i]) for i in range(self.num_generated)])
        self.predicted_ruwe = 10**pred_log_ruwe
        pred_sigma = np.concatenate([f_sigma(self.projected_sep[i], delta_g[i]) for i in range(self.num_generated)])

        # print(f"pred_log_ruwe: {np.shape(pred_log_ruwe)}")
        # print(f"pred_sigma: {np.shape(pred_sigma)}")
        # print(f"log_ruwe: {np.shape(log_ruwe)}")

        # g. Determine rejection probabilities
       
        # TODO: Not sure if I should change this either
        # rejection_prob = [stats.halfnorm.cdf(10**pred_log_ruwe[i], loc=10**log_ruwe, scale=10**pred_sigma[i]) for i in range(self.num_generated)]
        rejection_prob = stats.halfnorm.cdf(10**pred_log_ruwe, loc=10**log_ruwe, scale=10**pred_sigma)


        # never reject something where the observed ruwe is higher than the predicted (halfnorm)
        # never reject something that is outside the RUWE Distribution grid
       
        # TODO: Can be made into an array using bitwise comparisons (like we've done before)
        # rejection_prob = [0.0 if (not -.1 < delta_g[i] < 7.1) or
        #                  (not np.min(self.ruwe_dist['Sep(AU)']) < self.projected_sep[i] < np.max(self.ruwe_dist['Sep(AU)']))
        #                  else rejection_prob[i] for i in range(self.num_generated)]
        
        rejection_prob[(delta_g < -0.1) | (delta_g > 7.1)] = 0        
        
        rejection_prob[(self.projected_sep < np.min(self.ruwe_dist['Sep(AU)'])) | (self.projected_sep > np.max(self.ruwe_dist['Sep(AU)']))] = 0
        
        self.rejection_prob = rejection_prob

        rejection = np.random.rand(self.num_generated)
        
        # print(f"shape of rejection_prob: {np.shape(rejection_prob)}")
        # print(f"shape of rejection: {np.shape(rejection)}")

        
        # reject_list = [True if rejection[i] < rejection_prob[i] else False for i in range(self.num_generated)]
        # If we know how long something will be, a numpy array is faster.
        # Otherwise, regular lists are faster
        reject_list = rejection < rejection_prob



        return reject_list

    def load_stellar_model(self, star_age):
        # Read in file containing stellar model
        # TODO: Interpolate to get the chart for the exact age or binned age or something that doesnt rely on the age being in the chart
        model_chart = {}
        BHAC_file = f'{os.path.join(repo_path, "code/BHAC15_CFHT.txt")}'
        with open(BHAC_file, 'r') as content_file:
            content = content_file.read()
        tables = content.split(
            sep='\n----------------------------------------------------------------------------------------------------------------------------------------\n')
        tables = [x for x in tables if len(x) > 1]

        for table in tables:
            n = table.find('M')
            time_segment = table[0:n]
            table_segment = table[n:]
            age = float(time_segment[time_segment.find('=') + 1:])
            year_chart = Table.read(table_segment, format='ascii', fast_reader=False)
            year_chart['Age(Gyr)'] = age
            model_chart[age] = year_chart

        # Set up interpolation
        ages = np.array(list(model_chart.keys()))
        #  check if age is modeled
        if star_age in ages:
            return model_chart[star_age]
        #  find ages above and below the desired age
        diff = [star_age-x for x in ages]
        low_age = ages[np.argmin([x if x > 0 else float('inf') for x in diff])]
        high_age = ages[np.argmin([abs(x) if x < 0 else float('inf') for x in diff])]
        #  get masses
        new_model = Table()
        new_model['M/Ms'] = model_chart[high_age]['M/Ms']
        for y in model_chart[low_age].colnames[1:]:
            y_list = []
            for m in list(new_model['M/Ms']):
                # interpolate
                low_i = list(model_chart[low_age]['M/Ms']).index(m)
                high_i = list(model_chart[high_age]['M/Ms']).index(m)
                xs = [low_age, high_age]
                ys = [model_chart[low_age][y][low_i], model_chart[high_age][y][high_i]]
                f = scipy.interpolate.interp1d(xs, ys, kind='linear')
                y_list.append(f(star_age))
            #  add to table
            new_model[y] = y_list

        return new_model

    def get_gaia_info(self):
        if np.isfinite(self.gmag):
            #    Check if magnitude, color and Gaia solution are valid for calculating RUWE
            if 3.6 <= self.gmag <= 21. and -1 <= self.color <= 10:
                # All is well
                return
            else:
                # The magnitude or color is outside of bounds
                return -53

        else:
            #   Query Gaia for parallax, g_mag, color, astrometric_chi2 and n_good_obs, calculate RUWE
            coordinate = SkyCoord(self.star_ra, self.star_dec, frame='icrs')
            width = u.Quantity(10, u.arcsecond)
            height = u.Quantity(10, u.arcsecond)
            job_str = ("SELECT TOP 10 DISTANCE(POINT('ICRS', %f, %f), POINT('ICRS', ra, dec)) AS dist, * FROM gaiaedr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', %f, %f),CIRCLE('ICRS', ra, dec, 0.08333333)) ORDER BY dist ASC)" % (
                coordinate.ra.degree, coordinate.dec.degree, coordinate.ra.degree, coordinate.dec.degree))
            job = Gaia.launch_job(job_str)
            gaia_info = job.get_results()
            logging.info(f"gaia_info type: {type(gaia_info)}")
            # gaia_info.show_in_browser(jsviewer=True)
            if gaia_info and len(gaia_info) >= 1:
                # Only one star fits the coordinates, all is well
                self.gmag = float(gaia_info['phot_g_mean_mag'][0])
                self.gaia_id = int(gaia_info['source_id'][0])
                self.color = float(gaia_info['bp_rp'][0])
                self.n_good_obs = float(gaia_info['astrometric_n_good_obs_al'][0])
                self.astrometric_chi2 = float(gaia_info['astrometric_chi2_al'][0])
                self.parallax = float(gaia_info['parallax'][0])
                self.parallax_error = float(gaia_info['parallax_error'][0])
                self.ln_ruwe = float(np.log(gaia_info['ruwe'][0]))
                logging.info(f"gaia information:\n{self.gmag}\n{self.color}\n{self.n_good_obs}\n{self.astrometric_chi2}\n{self.parallax}\n{self.parallax_error}\n{self.ln_ruwe}\n")

                #    Check if magnitude, color and Gaia solution are valid for calculating RUWE
                if 3.6 <= self.gmag <= 21. and -1 <= self.color <= 10 and gaia_info['astrometric_params_solved'][0] == 31:
                    # All is well
                    return
                else:
                    # The magnitude or color is outside of bounds
                    return -53
            else:
                # No Gaia results
                logging.info("\nNo Gaia results!!!\n")
                return -51

    def read_dist(self):
        file_name = f'{os.path.join(repo_path, "code/RuweTableGP.txt")}'
        t = Table.read(file_name, format='ascii', delimiter=' ')

        self.ruwe_dist = t
        return
