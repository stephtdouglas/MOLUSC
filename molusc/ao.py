from datetime import datetime as dt
import datetime
import warnings
import os, pathlib, sys
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

warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.LinAlgWarning)

today = dt.today().isoformat().split("T")[0]

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
from molusc.utils import get_pro_sep, calc_anomaly_gaia

class AO:

    # Input file must be in mas and mags
    def __init__(self, filename, companions, star_mass, star_age, 
                 filter, obs_date=None, gaia=False):
        self.ao_filename = filename
        self.mass_ratio = companions.mass_ratio
        self.a = companions.a
        self.companions = companions
        self.star_mass = star_mass
        # Get the appropriate model and year
        self.age_model = self.load_stellar_model(filter, star_age)
        self.obs_date = obs_date
        if gaia:
            self.a_type = 'gaia'

    def analyze(self):
        num_generated = len(self.mass_ratio)

        # Unpack companions' orbital parameters
        period = self.companions.P
        e = self.companions.ecc
        arg_peri = self.companions.arg_peri
        phase = self.companions.phase
        cos_i = self.companions.cos_i
        if self.obs_date:
            T_0 = np.full(num_generated,fill_value=self.obs_date)
        else:
            T_0 = np.full(num_generated,fill_value=2457388.5)  # epoch 2016.0 in JD

        # Read in the contrast
        contrast = self.contrast
        a_type = self.a_type

        #todo add error handling for cases of too young and too old

        # Use interior function to calculate the model contrast
        self._calc_model_contrast()

        all_stars = np.stack([T_0,period,phase,e,arg_peri,cos_i,self.a],axis=1)
        
        try:
            cpu_ct = len(os.sched_getaffinity(0))
            print(f"Current time: {datetime.datetime.now()} -- AO cpu_count HPC:", cpu_ct)
        except AttributeError:
            cpu_ct = mp.cpu_count()-1
            print(f"Current time: {datetime.datetime.now()} -- AO cpu_count PC:", cpu_ct)
            
        divisor = num_generated // cpu_ct
        
        with Pool(cpu_ct) as pool:
            pro_sep = pool.starmap(get_pro_sep, all_stars, chunksize=divisor)

        four_arc = round(self.star_distance * 0.0000193906, 1)  # 4" in AU at distance of primary
        
        if a_type == 'hard limit':
            # Find Delta Mag limits, and/or recovery fraction
            # Interpolate linearly between given points to get the estimated contrast limit
            # Reject or accept the hypothetical binary based on the "hard limit" of the experimental contrast
            # 100% of binaries with a contrast less than the experimental contrast are rejected, 100% of those with
            # a greater contrast cannot be rejected
            
            contrast.sort('Sep (AU)')

            f_con = scipy.interpolate.interp1d(contrast['Sep (AU)'], contrast['Contrast'], kind='linear', bounds_error=False, fill_value=0)
                    
            # # Parallelziation
            # with Pool(cpu_ct) as pool:
            #     contrast_limit = pool.map(f_con, pro_sep, chunksize=divisor)

            contrast_limit = f_con(pro_sep)
                
                
            # End parallelization
                
            # Compare the model_contrast to the experimental_delta_K
            # If the model contrast is less than the experimental contrast it would have been seen and can be rejected
            # model contrast < contrast limit = reject (reject list = true)
            self.reject_list = np.greater(contrast_limit, self.model_contrast)#, dtype=bool)
            logging.info(self.reject_list)
            
        elif a_type == 'gradient':
            # Commented out just to be sure it doesn't go here for now :)
            pass
            # # Find Delta Mag limits, and/or recovery fraction
            # # Interpolate linearly between given points to get the estimated contrast limit
            # # Reject or accept the hypothetical binary based on a gradient of recovery rates for separation and contrast
            # recovery_rate = [0.]*num_generated

            # # Get column names and recovery rates
            # column_rates = [float(x.strip('%')) / 100.0 for x in list(contrast.columns)[1:]]
            # column_names = contrast.colnames[1:]

            # f = scipy.interpolate.interp2d(contrast['Sep (AU)'], column_rates, [contrast[x] for x in column_names])

            # for i in range(num_generated):
            #     column_names = contrast.colnames  # reset the list of column names
            #     # Interpolate
            #     if pro_sep[i] < contrast['Sep (AU)'][0]:  # closer than lowest limit, recovery rate = 0
            #         recovery_rate[i] = 0.
            #         continue
            #     elif pro_sep[i] > contrast['Sep (AU)'][-1]:  # further than farthest limit, recovery rate = 1
            #         recovery_rate[i] = 1.
            #         continue
            #     else:
            #         new_row = Table(rows=[[float(f(pro_sep[i], x)) for x in column_rates]], names=column_names[1:])
            #         new_row['Sep (AU)'] = pro_sep[i]
            #         # Determine which recovery rates the magnitude falls between, and assign it the lower one
            #         if self.model_contrast[i] < new_row[column_names[-1]]:  # Greater than largest recovery rate
            #             recovery_rate[i] = column_rates[-1]
            #             continue
            #         elif self.model_contrast[i] > new_row[column_names[1]]:  # Less than smallest recovery rate
            #             recovery_rate[i] = 0.
            #             continue
            #         else:
            #             for j in range(1, len(column_names)):
            #                 if new_row[column_names[j]][0] < self.model_contrast[i]:
            #                     recovery_rate[i] = column_rates[j-2]
            #                     break

            # # Make Reject list
            # random = np.random.uniform(0, 1, num_generated)
            # self.reject_list = [True if random[i] < recovery_rate[i] else False for i in range(0, num_generated)]

        # Write out information for display or output files
        self.pro_sep = pro_sep

        return np.array(self.reject_list)

    def analyze_gaia(self, gaia_limit):
        # Unpack companions' orbital parameters
        period = self.companions.P
        ecc = self.companions.ecc
        arg_peri = self.companions.arg_peri
        phase = self.companions.phase
        cos_i = self.companions.cos_i
        num_generated = len(self.mass_ratio)

        # Set date
        T_0 = 2457388.5  # epoch 2016.0 in JD

        # Read in the contrast
        contrast = self.contrast
        a_type = self.a_type

        # Find Delta Mag
        self._calc_model_contrast()

        # Calculate projected separation for each generated companion

        pro_sep = calc_anomaly_gaia(num_generated,T_0,period,phase,ecc,arg_peri,cos_i,self.a)
        
        #  Determine Gaia completeness detection limits
        four_arc = round(self.star_distance * 0.0000193906, 1)  # 4" in AU at distance of primary
        completness_absolute = gaia_limit - 5 * np.log10(self.star_distance / 2062650)  # apparent converted to absolute
        completeness_mag = np.round(completness_absolute - self.star_model_mag, 2)  # delta mag between the primary and gaia's detection limit

        # Adjusting to the gaia completeness mag
        # This seems to make just one value for completeness_mag, 
        # so no need to do a loop
        for colname in contrast.dtype.names[1:]:
            too_faint = contrast[colname] > completeness_mag
            contrast[colname][too_faint] = completeness_mag

        # Comment this section if you want nearest neighbor limits
        contrast['100%'] = np.zeros(len(contrast))
        contrast['100%'][contrast['Sep (AU)']>=four_arc] = completeness_mag
        new_row = np.full(len(contrast.colnames),fill_value=completeness_mag)
        new_row[0] = four_arc
        contrast.add_row(new_row)
        new_row[0] = self.star_distance
        contrast.add_row(new_row)
        contrast.sort('Sep (AU)')

        #  end neighbor-less segment

        # Uncomment this section if you want nearest neighbor limits
        # if self.nearest_neighbor_dist < four_arc:
        #     # If the nearest neighbor is less than 4 arcseconds away, I need to not add a 4" row, and truncate
        #     # the existing rows to a maximum of the nearest neighbor distance
        #     # first, find out what the interpolated limits are at the distance of the nearest neighbor
        #     column_rates = [float(x.strip('%')) / 100.0 for x in list(contrast.columns)[1:]]
        #     column_names = contrast.colnames[1:]
        #
        #     f_neighbor = scipy.interpolate.interp2d(contrast['Sep (AU)'], column_rates, [contrast[x] for x in column_names])
        #
        #     l = [round(float(f_neighbor(self.nearest_neighbor_dist, x)),2) for x in column_rates]
        #     contrast.add_row(([self.nearest_neighbor_dist]+l))
        #
        #     # Sort by separation
        #     contrast.sort('Sep (AU)')
        #     # Remove all rows after the nearest neighbor. There probably won't be any, since I'm assuming no
        #     # one is going to give me a contrast curve that has another star in it, but the code is here anyways
        #     ind = list(contrast['Sep (AU)']).index(self.nearest_neighbor_dist)
        #     contrast.remove_rows(slice(ind+1, len(contrast)))
        #
        # else:
        #     # Add two rows at the bottom, reaching from 4" to the nearest neighbor
        #     contrast['100%'] = [0.0]*len(contrast)
        #     contrast.add_row(([four_arc]+[completeness_mag]*(len(contrast.colnames)-1)))
        #     contrast.add_row(([self.nearest_neighbor_dist]+[completeness_mag]*(len(contrast.colnames)-1)))
        # # end nearest neighbor segment

        # Get column names and recovery rates
        # Column headers are recovery rates as percentages
        column_rates = [float(x.strip('%')) / 100.0 for x in list(contrast.columns)[1:]]
        column_names = contrast.colnames[1:]
        recovery_rate = np.zeros(num_generated)

        # TODO: interp2d is deprecated
        # 2D interpolation for contrast rate as a function of
        # separation and recovery rate
        #print(np.shape(contrast),len(contrast.columns))
        contr_map = np.asarray(contrast[column_names])
        contr_map = np.full((len(contrast),len(column_names)),fill_value=99.9)
        for i,colname in enumerate(column_names):
            contr_map[:,i] = contrast[colname]
        #contr_map = np.vstack([np.asarray(row) for row in contrast[column_names]])
        calc_contr = scipy.interpolate.RegularGridInterpolator(
            (np.asarray(contrast['Sep (AU)']),column_rates),
            contr_map,
            bounds_error=False,fill_value=np.inf)

        # closer than lowest limit, recovery rate = 0
        recovery_rate[pro_sep<contrast['Sep (AU)'][0]] = 0

        # further than farthest limit, recovery rate = 1
        recovery_rate[pro_sep>contrast['Sep (AU)'][-1]] = 1

        intermediate_contrast = np.where((pro_sep>=contrast['Sep (AU)'][0]) |
                                        (pro_sep<=contrast['Sep (AU)'][-1]))[0]

        # Within the appropriate limits, interpolate to find
        # the recovery rate
        nrates = len(column_rates)
        for i in intermediate_contrast:
#            column_names = contrast.colnames  # reset the list of column names
            # Interpolate
            pair_contr = np.full((nrates,2),pro_sep[i])
            pair_contr[:,1] = column_rates
            rec_contr = calc_contr(pair_contr)

            # Determine which recovery rates the magnitude falls between, and assign it the lower one
            j = np.where(rec_contr<=self.model_contrast[i])[0]
            if len(j)==0:
                # the model contrast is fainter than any recovery rate limit
                recovery_rate[i] = 0
            else:
                recovery_rate[i] = column_rates[j[0]]

        # Make Reject list
        # TODO: why is this rejection selection different from rv?
        random = np.random.uniform(0, 1, num_generated)
        self.reject_list = random < recovery_rate

        # Write out information for display or output files
        self.pro_sep = pro_sep

        return np.array(self.reject_list)

    def _calc_model_contrast(self):

        # Determine low and high mass limits
        self.low_mass_limit = self.age_model['M/Ms'][0]
        self.high_mass_limit = self.age_model['M/Ms'][-1]

        if self.star_mass > self.high_mass_limit:
            return -23
        elif self.star_mass < self.low_mass_limit:
            return -24

        # Find model mag of primary star
        self.star_model_mag = self.find_mag(self.star_mass, self.age_model)
        #self.star_model_mag = star_model_mag
        print(f'Current time: {datetime.datetime.now()} -- Star Model Mag {self.star_model_mag}')

        # Get masses of companion stars
        # TODO: why are we rounding this?
        cmp_mass = np.round(self.star_mass * self.mass_ratio , 3) # companion mass in solar masses
        # del(self.mass_ratio)

        # Get companion star magnitudes, assign infinite magnitude if below lowest modeled mass
        f_mag = scipy.interpolate.interp1d(self.age_model['M/Ms'], self.age_model['Mag'], kind='cubic', 
                                           fill_value=np.nan, bounds_error=False)
        cmp_model_mag = f_mag(cmp_mass)

        # Find Delta Mag
        self.model_contrast = cmp_model_mag - self.star_model_mag



    def find_mag(self, star_mass, t):
        f_mag = scipy.interpolate.interp1d(t['M/Ms'], t['Mag'], kind='cubic', fill_value='extrapolate')
        mag = f_mag(star_mass)
        return mag


    def get_distance(self, star_RA, star_DEC, parallax=np.nan):
        # coordinate = SkyCoord(star_RA, star_DEC, frame='icrs')
        # width = u.Quantity(10, u.arcsecond)
        # height = u.Quantity(10, u.arcsecond)
        # job_str = ("SELECT TOP 10 DISTANCE(POINT('ICRS', %f, %f), POINT('ICRS', ra, dec)) AS dist, * FROM gaiaedr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', %f, %f),CIRCLE('ICRS', ra, dec, 0.08333333)) ORDER BY dist ASC)" % (coordinate.ra.degree, coordinate.dec.degree, coordinate.ra.degree, coordinate.dec.degree))
        # job = Gaia.launch_job(job_str)
        # gaia_info = job.get_results()
        # gaia_info.show_in_browser(jsviewer=True)
        # print(f"----------------------\nHere is the gaia info table thing: {gaia_info}\n----------------------")

        if np.isfinite(parallax):
            # Convert star parallax to distance in AU: d[AU] = 1/p["] * 206265 AU/parsec
            star_distance = 1 / (parallax / 1000) * 206265
            self.star_distance = star_distance
            self.nearest_neighbor_dist = np.nan

        else:
            coordinate = SkyCoord(star_RA, star_DEC, frame='icrs')
            width = u.Quantity(10, u.arcsecond)
            height = u.Quantity(10, u.arcsecond)
            job_str = ("SELECT TOP 10 DISTANCE(POINT('ICRS', %f, %f), POINT('ICRS', ra, dec)) AS dist, * FROM gaiaedr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', %f, %f),CIRCLE('ICRS', ra, dec, 0.08333333)) ORDER BY dist ASC)" % (coordinate.ra.degree, coordinate.dec.degree, coordinate.ra.degree, coordinate.dec.degree))
            job = Gaia.launch_job(job_str)
            gaia_info = job.get_results()

            if gaia_info: # Something in here has the gaia id! maybe called "source id". Find this and save as self.gaia_id
                if len(gaia_info) > 1:
                    # Multiple possible sources. Sort by search distance and take the closest one. Print warning.
                    # Nearest neighbor distance is set to distance between chosen source and it's closest neighbor in the search region
                    gaia_info.sort('dist')
                    parallax = gaia_info['parallax'][0]
                    picked_coord = SkyCoord(gaia_info['ra'][0], gaia_info['dec'][0], unit='degree')
                    # Convert star parallax to distance in AU: d[AU] = 1/p["] * 206265 AU/parsec
                    star_distance = 1 / (parallax / 1000) * 206265
                    self.star_distance = star_distance
                    self.gaia_id = gaia_info["source_id"][0]

                    # Find nearest neighbor
                    #    calculate distance between the picked star and all other stars in search radius
                    gaia_info['separation'] = [picked_coord.separation(SkyCoord(x['ra'], x['dec'], unit='degree')).mas for x in gaia_info]
                    #    sort the search results by separation from picked star (the picked star should have a separation of zero)
                    gaia_info.sort('separation')
                    #    choose the closest one to the picked star
                    nearest_neighbor_dist = gaia_info['separation'][1]  # distance to n.n. in mas
                    #   convert from mas to AU
                    self.nearest_neighbor_dist = round(self.star_distance * np.tan(np.radians(nearest_neighbor_dist/(3.6e6))), 1)  # distance to n.n. in AU
                    return -52
                else:
                    # Only once source in search area. Nearest neighbor distance is set to search width
                    parallax = gaia_info['parallax'][0]
                    # Convert star parallax to distance in AU: d[AU] = 1/p[as] * 206265AU/parsec
                    star_distance = 1 / (parallax / 1000) * 206265
                    self.star_distance = star_distance
                    self.gaia_id = gaia_info["source_id"][0]


                    # Set nearest neighbor distance to search distance
                    nearest_neighbor_dist = width.to('mas').value
                    #   convert from mas to AU
                    self.nearest_neighbor_dist = round(self.star_distance * np.tan(np.radians(nearest_neighbor_dist/(3.6e6))), 1)  # distance to n.n. in AU
                    return 0
            else:
                return -51

    def load_stellar_model(self, filter, star_age):
        print('LOADING STELLAR MODEL...')
        # Read in file containing stellar model with the filter needed
        if filter == 'J' or filter == 'H' or filter == 'K':  # 2MASS filters
            model_chart = {}
            BHAC_file = f'{os.path.join(repo_path, "reference_data/BHAC15_2MASS.txt")}'
            with open(BHAC_file, 'r') as content_file:
                content = content_file.read()
            tables = content.split(
                sep='\n-----------------------------------------------------------------------------------------------\n')
            tables = [x for x in tables if len(x) > 1]

            for table in tables:
                n = table.find('M')
                time_segment = table[0:n]
                table_segment = table[n:]
                age = float(time_segment[time_segment.find('=') + 1:])
                year_chart = Table.read(table_segment, format='ascii', fast_reader=False, comment='!')
                if filter == 'J':
                    year_chart = year_chart['M/Ms', 'Mj']
                    year_chart.rename_column('Mj', 'Mag')
                elif filter == 'H':
                    year_chart = year_chart['M/Ms', 'Mh']
                    year_chart.rename_column('Mh', 'Mag')
                elif filter == 'K':
                    year_chart = year_chart['M/Ms', 'Mk']
                    year_chart.rename_column('Mk', 'Mag')
                model_chart[age] = year_chart
        elif filter == 'R' or filter == 'I':
            model_chart = {}
            BHAC_file = f'{os.path.join(repo_path, "reference_data/BHAC15_CFHT.txt")}'
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
                year_chart = Table.read(table_segment, format='ascii', fast_reader=False, comment='!')
                if filter == 'R':
                    year_chart = year_chart['M/Ms', 'R']
                    year_chart.rename_column('R', 'Mag')
                elif filter == 'I':
                    year_chart = year_chart['M/Ms', 'I']
                    year_chart.rename_column('I', 'Mag')
                model_chart[age] = year_chart
        elif filter == 'G' or filter == 'Rp' or filter ==  'Bp':
            model_chart = {}
            BHAC_file = f'{os.path.join(repo_path, "reference_data/BHAC15_GAIA.txt")}'
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
                year_chart = Table.read(table_segment, format='ascii', fast_reader=False, comment='!')
                if filter == 'G':
                    year_chart = year_chart['M/Ms', 'G']
                    year_chart.rename_column('G', 'Mag')
                elif filter == 'Rp':
                    year_chart = year_chart['M/Ms', 'G_RP']
                    year_chart.rename_column('G_RP', 'Mag')
                elif filter == 'Bp':
                    year_chart = year_chart['M/Ms', 'G_BP']
                    year_chart.rename_column('G_BP', 'Mag')
                model_chart[age] = year_chart
        elif filter == 'L' or filter == 'LL' or filter == 'M':
            print('MODEL:  CIT2')
            model_chart = {}
            BHAC_file = f'{os.path.join(repo_path, "reference_data/BHAC15_CIT2.txt")}'
            with open(BHAC_file, 'r') as content_file:
                content = content_file.read()
            content = content[content.find('\n\n\n'):]  # Cut off the intro material
            tables = content.split(
                sep='\n!-----------------------------------------------------------------------------------------------\n\n')
            tables = [x for x in tables if len(x) > 1]
            for table in tables:
                n1 = table.find('t (Gyr)')
                n2 = table.find('M')
                time_segment = table[n1:table[n1:].find('!') + n1]
                table_segment = table[n2:]

                age = float(time_segment[time_segment.find('=') + 1:])
                year_chart = Table.read(table_segment, format='ascii', fast_reader=False, comment='!')

                if filter == 'L':
                    year_chart = year_chart['M/Ms', 'Ml']
                    year_chart.rename_column('Ml', 'Mag')
                elif filter == 'LL':
                    year_chart = year_chart['M/Ms', 'Mll']
                    year_chart.rename_column('Mll', 'Mag')
                elif filter == 'M':
                    year_chart = year_chart['M/Ms', 'Mm']
                    year_chart.rename_column('Mm', 'Mag')
                model_chart[age] = year_chart
        # Set up interpolation
        ages = np.array(list(model_chart.keys()))
        #  check if age is modeled, if it is simply take that chart, if not interpolate within the model to get it
        if star_age in ages:
            return model_chart[star_age]
        #  find ages above and below the desired age
        diff = [star_age - x for x in ages]
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
        del(y_list)
        return new_model

    def read_contrast(self):
        contrast = {}
        # Each line in the file should have the separation in mas and the contrast for different recovery rates,
        # separated by whitespace
        # If a date is included it should be the first line of the file, starting with a # and followed by the JD date
        if os.path.exists(self.ao_filename)==False:
            return -21

        try:
            with open(self.ao_filename, 'r') as f:
                read_data = f.read()
                if read_data.startswith('#'):
                    # Remove the date from the first part
                    split = read_data.find('\n')
                    # Get the date
                    self.obs_date = float(read_data[1:split].lstrip())
                    # Get the contrast
                    contrast = read_data[split:]
                    contrast_table = Table.read(contrast, format='ascii')
                else:
                    # No date given, just get the contrast
                    contrast_table = Table.read(read_data, format='ascii', delimiter=' ', fast_reader=False)
                # Rename columns, assuming they are in the correct format of separation, magnitude
                contrast_table.rename_column(list(contrast_table.columns)[0], 'Sep')
                # Convert separation from mas to AU, order columns correctly
                contrast_table['Sep (AU)'] = [round(self.star_distance * np.tan(np.radians(x/(3.6e6))), 1) for x in contrast_table['Sep']]
                order = ['Sep (AU)'] + list(contrast_table.columns)[1:-1]
                contrast_table = contrast_table[order]
        except TypeError:
            return -22

        # If there is only one contrast given (hard limit)
        if len(contrast_table.columns) == 2:
            self.a_type = 'hard limit'
            contrast_table.rename_column(list(contrast_table.columns)[1], 'Contrast')
        # Gradient contrast given
        # The first line should be a header, and should list the recovery rates used
        else:
            if self.a_type != 'gaia':
                self.a_type = 'gradient'

        self.contrast = contrast_table

        return 0



class AO_detection(AO):
    # Inheriting: get_distance, _calc_model_contrast
    # technically also read_contrast but there's no point in using it here

    def __init__(self, companions, star_mass, star_age, 
                 filter, pa, e_pa, rho, e_rho, dM, e_dM,
                 obs_date=None):
        """
        Given an observed position angle (pa), angular separation (rho),
        delta mag (dM), and observation filter, reject possible companions

        pa, rho, dM, and filter can be numbers or arrays (all must be the same
        type and length if appropriate). If multiple detections are provided,
        then each companion in the prior will be checked against all of them
        """
        self.pa = pa
        self.e_pa = e_pa
        self.rho = rho 
        self.e_rho = e_rho 
        self.dM = dM 
        self.e_dM = e_dM
        self.filter = filter
        self.companions = companions
        self.obs_date = obs_date
        self.star_mass = star_mass
        self.star_age = star_age

    def analyze(self):
        
        num_generated = comps.num_generated

        # Unpack companions' orbital parameters
        period = self.companions.P
        e = self.companions.ecc
        arg_peri = self.companions.arg_peri
        phase = self.companions.phase
        cos_i = self.companions.cos_i

        # TODO: we should always care about the observation date
        if self.obs_date:
            T_0 = np.full(num_generated,fill_value=self.obs_date)
        else:
            T_0 = np.full(num_generated,fill_value=2457388.5)  # epoch 2016.0 in JD

        # Get the appropriate age model(s)
        # TODO: right now this only assumes one filter at a time. 
        # I can either run on one filter at a time for consistency, or do them both together

        self._calc_model_contrast()

        # First rejection step - contrast of the prior companion doesn't match

        max_dM = dM + e_dM * 3
        min_dM = dM - e_dM * 3

        if isinstance(dM, (list, tuple, np.ndarray)):

