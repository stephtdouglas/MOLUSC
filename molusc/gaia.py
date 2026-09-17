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
from time import time
from astropy.utils.exceptions import AstropyWarning
from multiprocessing import Process
from multiprocessing.pool import Pool
import multiprocessing as mp
# from guppy import hpy
# import tracemalloc

# from pkgcore.config import load_config
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.LinAlgWarning)

# c = load_config()
# hp = hpy()
today = dt.today().isoformat().split("T")[0]

import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent
from molusc.utils import get_pro_sep, calc_anomaly_gaia
from molusc.ao import AO

class Gaia(AO):

    def analyze_gaia(self, gaia_limit,nearest_neighbor=False):
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

        # Determine low and high mass limits
        low_mass_limit = self.age_model['M/Ms'][0]
        high_mass_limit = self.age_model['M/Ms'][-1]

        if self.star_mass > high_mass_limit:
            return -23
        elif self.star_mass < low_mass_limit:
            return -24

        # Find model mag of primary star
        self.star_model_mag = self.find_mag(self.star_mass, self.age_model)
        print(f'Current time: {datetime.datetime.now()} -- Star Model Mag {self.star_model_mag}')  # TESTING

        # Get masses of companion stars
        cmp_mass = self.star_mass * self.mass_ratio  # companion mass in solar masses
        cmp_mass = np.round(cmp_mass, 3)

        # Get companion star magnitudes, assign infinite magnitude if below lowest modeled mass
        f_mag = scipy.interpolate.interp1d(self.age_model['M/Ms'], self.age_model['Mag'], 
                                           kind='cubic', fill_value=np.inf, bounds_error=False)
        cmp_model_mag = f_mag(cmp_mass)

        # Find Delta Mag
        model_contrast = cmp_model_mag - self.star_model_mag

        # Calculate projected separation for each generated companion

        pro_sep = calc_anomaly_gaia(num_generated,T_0,period,phase,ecc,arg_peri,cos_i,self.a)
        
        #  Determine Gaia completeness detection limits
        # at Praesepe's distance this gives 4" -> 10^7 AU, like 10 pc??
        four_arc = round(self.star_distance * 0.0000193906, 1)  # 4" in AU at distance of primary
        completness_absolute = gaia_limit - 5 * np.log10(self.star_distance / 2062650)  # apparent converted to absolute
        completeness_mag = np.round(completness_absolute - self.star_model_mag, 2)  # delta mag between the primary and gaia's detection limit

        # Adjusting to the gaia completeness mag
        # This seems to make just one value for completeness_mag, 
        # so no need to do a loop
        for colname in contrast.dtype.names[1:]:
            too_faint = contrast[colname] > completeness_mag
            contrast[colname][too_faint] = completeness_mag

        # Skip nearest neighbor analysis
        if nearest_neighbor is False:
            contrast['100%'] = np.zeros(len(contrast))
            contrast['100%'][contrast['Sep (AU)']>=four_arc] = completeness_mag
            new_row = np.full(len(contrast.colnames),fill_value=completeness_mag)
            new_row[0] = four_arc
            # This adds a new row at an insanely high separation
            # at Praesepe's distance it's 10^7 AU, like 10 pc??
            # But it does sort it into the right place
            contrast.add_row(new_row)
            new_row[0] = self.star_distance
            contrast.add_row(new_row)
            contrast.sort('Sep (AU)')
        # Run nearest neighbor analysis
        elif (nearest_neighbor is True) and (self.nearest_neighbor_dist <= four_arc):
            # If the nearest neighbor is less than 4 arcseconds away, I need to not add a 4" row, and truncate
            # the existing rows to a maximum of the nearest neighbor distance
            # first, find out what the interpolated limits are at the distance of the nearest neighbor
            column_rates = [float(x.strip('%')) / 100.0 for x in list(contrast.columns)[1:]]
            column_names = contrast.colnames[1:]
        
            contr_map = np.asarray(contrast[column_names])
            contr_map = np.full((len(contrast),len(column_names)),fill_value=99.9)
            for i,colname in enumerate(column_names):
                contr_map[:,i] = contrast[colname]
            f_neighbor = scipy.interpolate.RegularGridInterpolator((np.asarray(contrast['Sep (AU)']),column_rates),
                                                                    contr_map)#,bounds_error=False,fill_value=np.inf)
        
            l = [round(float(f_neighbor(self.nearest_neighbor_dist, x)),2) for x in column_rates]
            contrast.add_row(([self.nearest_neighbor_dist]+l))
        
            # Sort by separation
            contrast.sort('Sep (AU)')
            # Remove all rows after the nearest neighbor. 
            ind = list(contrast['Sep (AU)']).index(self.nearest_neighbor_dist)
            contrast.remove_rows(slice(ind+1, len(contrast)))
        
        elif (nearest_neighbor is True) and (self.nearest_neighbor_dist > four_arc):
            # Add two rows at the bottom, reaching from 4" to the nearest neighbor
            contrast['100%'] = [0.0]*len(contrast)
            contrast.add_row(([four_arc]+[completeness_mag]*(len(contrast.colnames)-1)))
            contrast.add_row(([self.nearest_neighbor_dist]+[completeness_mag]*(len(contrast.colnames)-1)))
        else:
            return -56
        # end nearest neighbor segment


        # Get column names and recovery rates
        # Column headers are recovery rates as percentages
        column_rates = [float(x.strip('%')) / 100.0 for x in list(contrast.columns)[1:]]
        column_names = contrast.colnames[1:]
        recovery_rate = np.zeros(num_generated)

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
        # because we're insensitive to it
        recovery_rate[pro_sep<contrast['Sep (AU)'][0]] = 0

        # further than farthest limit, recovery rate = 0
        # Because while it should be resolved in that case AND THIS IS GAIA
        recovery_rate[pro_sep>contrast['Sep (AU)'][-1]] = 1

        # Where we actually have limits, more work needed
        # This selection should be an and, not an or
        intermediate_contrast = np.where((pro_sep>=contrast['Sep (AU)'][0]) &
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
            j = np.where(rec_contr<=model_contrast[i])[0]
            if len(j)==0:
                # the model contrast is fainter than any recovery rate limit
                recovery_rate[i] = 0
            else:
                recovery_rate[i] = column_rates[j[0]]

        # Make Reject list
        # TODO: why is this rejection selection different from rv?
        random = np.random.uniform(0, 1, num_generated)
        # If recovery_rate is 1, the companion will always be rejected
        # If recovery_rate is 0, the companion will always be accepted
        self.reject_list = random < recovery_rate

        # Write out information for display or output files
        self.model_contrast = model_contrast
        self.pro_sep = pro_sep
        self.low_mass_limit = low_mass_limit

        return np.array(self.reject_list)


