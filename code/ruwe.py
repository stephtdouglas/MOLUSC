from datetime import datetime
import numpy as np
import scipy as scipy
import scipy.stats as stats
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
from astroquery.gaia import Gaia
import warnings
import os
from astropy.utils.exceptions import AstropyWarning
import logging
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

today = datetime.today().isoformat().split("T")[0]
global repo_path
repo_path = os.getenv('MOLOC').replace("\\", "/")

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
        self.num_generated = companions.get_num()
        self.a = companions.get_a()
        self.period = companions.get_P()
        self.e = companions.get_ecc()
        self.phase = companions.get_phase()
        self.arg_peri = companions.get_arg_peri()
        self.cos_i = companions.get_cos_i()
        self.mass_ratio = companions.get_mass_ratio()
        self.gmag = np.nan

    def analyze(self):
        # a. Calculate projected separation
        pro_sep = [0.0 for i in range(self.num_generated)]
        T_0 = 2457388.5  # epoch 2016.0 in JD, corresponding to Gaia eDR3

        for i in range(self.num_generated):
            # 1. Calculate mean anomaly
            M = 2 * np.pi * T_0 / self.period[i] - self.phase[i]
            # 2. Calculate eccentric anomaly iteratively
            prev_E = 0.0
            current_E = M
            while abs(current_E - prev_E) > 0.00001:
                prev_E = current_E
                current_E = M + self.e[i] * np.sin(prev_E)
            # 3. Calculate true anomaly
            f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - self.e[i]) / (1 + self.e[i])))
            # 4. Calculate projected separation in AU
            # Added an absolute value around the cos(f) in the separation calculation
            alpha = f + self.arg_peri[i]
            sqtmp = np.sqrt(np.sin(alpha)**2+np.cos(alpha)**2 * self.cos_i[i]**2)
            pro_sep[i] = self.a[i] * (1-self.e[i]**2)/(1+self.e[i]*np.cos(f))*sqtmp  # AU
        self.projected_sep = pro_sep

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

        pred_log_ruwe = np.concatenate([f_ruwe(self.projected_sep[i], delta_g[i]) for i in range(self.num_generated)])
        self.predicted_ruwe = 10**pred_log_ruwe
        pred_sigma = np.concatenate([f_sigma(self.projected_sep[i], delta_g[i]) for i in range(self.num_generated)])

        # g. Determine rejection probabilities
        rejection_prob = [stats.halfnorm.cdf(10**pred_log_ruwe[i], loc=10**log_ruwe, scale=10**pred_sigma[i]) for i in range(self.num_generated)]

        # never reject something where the observed ruwe is higher than the predicted (halfnorm)
        # never reject something that is outside the RUWE Distribution grid
        rejection_prob = [0.0 if (not -.1 < delta_g[i] < 7.1) or
                         (not np.min(self.ruwe_dist['Sep(AU)']) < self.projected_sep[i] < np.max(self.ruwe_dist['Sep(AU)']))
                         else rejection_prob[i] for i in range(self.num_generated)]
        self.rejection_prob = rejection_prob

        rejection = np.random.rand(self.num_generated)

        reject_list = [True if rejection[i] < rejection_prob[i] else False for i in range(self.num_generated)]

        return reject_list

    def load_stellar_model(self, star_age):
        # Read in file containing stellar model
        # todo Interpolate to get the chart for the exact age or binned age or something that doesnt rely on the age being in the chart
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
            gaia_info.show_in_browser(jsviewer=True)
            if gaia_info and len(gaia_info) >= 1:
                # Only one star fits the coordinates, all is well
                self.gmag = gaia_info['phot_g_mean_mag'][0]
                self.color = gaia_info['bp_rp'][0]
                self.n_good_obs = gaia_info['astrometric_n_good_obs_al'][0]
                self.astrometric_chi2 = gaia_info['astrometric_chi2_al'][0]
                self.parallax = gaia_info['parallax'][0]
                self.parallax_error = gaia_info['parallax_error'][0]
                self.ln_ruwe = np.log(gaia_info['ruwe'][0])
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
