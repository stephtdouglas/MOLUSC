from datetime import datetime as dt
import datetime
import warnings
import os, pathlib
import logging

import numpy as np
import scipy as scipy
import scipy.stats as stats
from astropy.table import Table
from time import time
import multiprocessing as mp
import gc
from astropy.utils.exceptions import AstropyWarning
from timeit import timeit
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

today = dt.today().isoformat().split("T")[0]
import molusc
repo_path = pathlib.Path(molusc.__file__).resolve().parent.parent

def calculate_RV_parallel(period, mass_ratio, a, ecc, cos_i, arg_peri, phase, MJD, calc):
    # Exactly the same as calculate_RV, but with an extra parameter stating whether you need to calculate RV
    # Calculates the RVs for each item when passed arrays of orbital parameters
    # Inputs: Arrays of Period, Mass Ratio, Semi-Major Axis, eccentricity, inclination, arg peri, phase, calculation times
    # Outputs: Velocity Semi-Amplitude (km/s), RVs at each time in MJD


    print(f'\n\nCurrent time: {datetime.datetime.now()} -- Calculating RVs Parallel function {mp.current_process()}...')
    sin_i = np.sin(np.arccos(cos_i))

    n = len(period)
    ndates = len(MJD)
    # Create a blank array with one row for every input period
    # and each row contains RVs equivalent to the observation dates
    RV = np.zeros((n,ndates))
    # RV = [[0.0 for i in range(len(MJD))] for j in range(n)]
    
    a_star = a * (mass_ratio / (mass_ratio + 1))
    # a_star = np.multiply(a, np.divide(mass_ratio, np.add(mass_ratio, 1)))

    K = (2*np.pi/period) * (a_star * sin_i) / np.sqrt(1-ecc**2)
    # K = np.multiply(np.divide((2 * np.pi), period),np.divide(np.multiply(a_star, sin_i), np.sqrt((1 - np.square(ecc)))))  # AU/days
    K = K * 1731.48  # km/s

    # TODO: Can this be converted to array operations for faster speed?
    for i in range(n):  # Iterate over companions
        if calc[i]:
            for j in range(0, len(MJD)):  # Iterate over times
                # Find E
                M = 2 * np.pi * MJD[j] / period[i] - phase[i]
                prev_E = 0.0
                current_E = M

                while abs(current_E - prev_E) > 0.00001:
                    prev_E = current_E
                    current_E = M + ecc[i] * np.sin(prev_E)

                # Find true anomaly, f
                f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - ecc[i]) / (1 + ecc[i])))
                # Find predicted RV
                RV[i][j] = K[i] * (np.sin(arg_peri[i] + f) + ecc[i] * np.sin(arg_peri[i]))  # km/s
    print(f'Current time: {datetime.datetime.now()} -- Finished calculating RVs! Parallel {mp.current_process()}')

    return K, RV

def zero_point_model(prediction, a):
    # Alters the prediction by some zero point shift, a
    y = prediction + a
    return y

def zero_point_fit_parallel(experimental, errors, predicted):

    A = np.zeros(len(predicted))
    for i in range(len(predicted)):
        [a], _ = scipy.optimize.curve_fit(zero_point_model, predicted[i], experimental, sigma=errors)
        A[i] = a
    return A

class RV:
    # Needs from application: filename, star mass, all generated values (P, e, arg_peri, cos_i, mass_ratio, a)
    # Needs to give back reject list
    # class variables
    rv_filename = ''
    added_jitter = 0.02 # km/s
    rv_floor = 0.02  # km/s
    rv_reject_list = []
    jitter_reject_list = []
    MJD = []
    experimental_RV = []
    measurement_error = []
    predicted_RV = []
    jitter = []

    # class functions
    def __init__(self, filename, resolution, companions, mass, age, added_jitter=0, rv_floor=20, extra_output=True):
        self.restore_defaults()
        self.rv_filename = filename
        self.resolution = resolution
        self.companions = companions
        self.star_mass = mass
        self.added_jitter = added_jitter/1000  # convert to km/s
        self.rv_floor = rv_floor/1000  # km/s
        self.extra_output = extra_output
        
        print(f'Current time: {datetime.datetime.now()} -- \n\nLoading stellar model...\n\n')
        self.model = self.load_stellar_model('G', age)
        print(f'Current time: {datetime.datetime.now()} -- Finished loading stellar model')

    def analyze_rv(self):
        # Calculate predicted RV
        #     Using orbital mechanics equations from Perryman 2011, and solving for E(anomaly) numerically
        #     For time comparison I want to generate times with the same range as the times in MJD, and then calculate
        #     the predicted RV at that time, which I can then compare to the experimental values using least squares
        num_generated = self.companions.num_generated
        # Unpack parameters
        print(f'Current time: {datetime.datetime.now()} -- Unpacking parameters...')
        period = self.companions.P
        mass_ratio = self.companions.mass_ratio
        a = self.companions.a
        ecc = self.companions.ecc
        cos_i = self.companions.cos_i
        arg_peri = self.companions.arg_peri
        phase = self.companions.phase
        print(f'Current time: {datetime.datetime.now()} -- Finished unpacking parameters')
        nobs = len(self.MJD)

        # Determine velocity limit
        print(f'Current time: {datetime.datetime.now()} -- Determining velocity limit...')
        delta_v = 2.998e5 / self.resolution  # m/s

        t = time()

        # Determine contrast
        cmp_mass = np.multiply(self.star_mass, mass_ratio)  # companion mass in solar masses

        f_mag = scipy.interpolate.interp1d(self.model['M/Ms'], self.model['Mag'], kind='cubic', fill_value='extrapolate')
        # self.model.show_in_browser(jsviewer=True, tableid="self.model")


        prim_model_mag = f_mag(self.star_mass)
        # companion mags, assign inf if below lowest modeled mass
        cmp_model_mag = np.full(num_generated,np.inf)
        apply_mag = cmp_mass>=self.model["M/Ms"][0]
        cmp_model_mag[apply_mag] = f_mag(cmp_mass[apply_mag])
        contrast = np.subtract(cmp_model_mag, prim_model_mag)
        print(f'Current time: {datetime.datetime.now()} -- Finished determining velocity limit')


        # Determine luminosity
        print(f'Current time: {datetime.datetime.now()} -- Determining luminosity...')
        f_lum = scipy.interpolate.interp1d(self.model['M/Ms'], self.model['L/Ls'], kind='cubic', fill_value='extrapolate')
        prim_lum = np.power(10, f_lum(self.star_mass))  # primary luminosity
        # companion luminosity, 0 if below limit
        # (same limit as for mag)
        cmp_lum = np.zeros(num_generated)
        cmp_lum[apply_mag] = np.power(10, f_lum(cmp_mass[apply_mag]))
        print(f'Current time: {datetime.datetime.now()} -- Determined luminosity')

        # Choose to run either parallelized or non-parallelized based on the number of generated companions
        print(f'Current time: {datetime.datetime.now()} -- Determining whether or not to run parallelized...')
        if num_generated < 50000:  # Serial
            print(f'Current time: {datetime.datetime.now()} -- Decided to run non-parallelzied')
            self.predicted_RV = [np.zeros(len(self.MJD)) for x in range(num_generated)] # pre-allocate

            # calculate RV curves
            prim_K, prim_rv = self.calculate_RV(period, mass_ratio, a, ecc, cos_i, arg_peri, phase, self.MJD)
            # cmp_K, cmp_rv = self.calculate_RV(period, np.divide(1, mass_ratio), a, ecc, cos_i, arg_peri, phase, self.MJD)
            cmp_K = prim_K / mass_ratio
            cmp_rv = -1 * prim_rv / mass_ratio
            # cmp_rv = np.multiply(-1, cmp_rv)

            max_delta_rv = np.max(np.absolute(np.subtract(prim_rv, cmp_rv)), axis=1)

            # Determine the overall predicted RV
            print(f'Current time: {datetime.datetime.now()} -- Determining overall predicted RV...')

            # TODO: this has probably become wrong somewhere along the way
            comb_rvs = np.hstack(prim_rv,cmp_rv)
            prim_lum_arr = np.full_like(cmp_lum,fill_value=prim_lum)
            comb_weights = np.hstack(prim_lum_arr, cmp_lum)
            self.predicted_RV = np.average(comb_rvs,axis=1,weights=comb_weights)
            # SB1s have high contrast
            sb1 = contrast > 5
            self.predicted_RV[sb1] = prim_rv[sb1]
            sb2 = np.where(contrast>=5)[0]

            
            print(f'Current time: {datetime.datetime.now()} -- Running zero point models...')
            # TODO: can this be further streamlined?
            for i in range(0, num_generated):
                # Fit the zero point
                [zero_point], pcov = scipy.optimize.curve_fit(self.zero_point_model, self.predicted_RV[i], self.experimental_RV, sigma=self.measurement_error)
                # Shift all predicted values by the zero_point
                self.predicted_RV[i] += zero_point
            print(f'Current time: {datetime.datetime.now()} -- Finished running zero point models')
            
            print(f'Current time: {datetime.datetime.now()} -- Determined overall predicted RV')

            print(f'Current time: {datetime.datetime.now()} -- Comparing experimental RV to predicted RV...')
            # Compare experimental and predicted RVs
            amp = np.ptp(self.predicted_RV, axis=1)
            # amp = [np.ptp(self.predicted_RV[i]) for i in range(num_generated)]

            model_difference = self.experimental_RV - self.predicted_RV
            sq_term = self.measurement_error**2 + self.added_jitter**2
            chi_squared0 = model_difference / sq_term
            chi_squared = np.sum(chi_squared0,axis=1)

            # The degrees of freedom is equal to (N-1)+1, for the number of data points and the applied velocity shift
            ndatesm1 = len(self.MJD)-1
            # TODO: can this be further streamlined?
            prob = [stats.chi2.cdf(chi_squared[i], ndatesm1) for i in range(num_generated)]
            print(f'Current time: {datetime.datetime.now()} -- Finished comparing experimental RV to predicted RV')

        else:  
            # Parallel
            print(f'Current time: {datetime.datetime.now()} -- Decided to run parallelzied')
            # Determine cpu count
            try:
                cpu_ct = len(os.sched_getaffinity(0))
                print(f"Current time: {datetime.datetime.now()} -- RV cpu_count HPC:", cpu_ct)
            except AttributeError:
                cpu_ct = mp.cpu_count()-1
                print(f"Current time: {datetime.datetime.now()} -- RV cpu_count PC:", cpu_ct)


            print(f'Current time: {datetime.datetime.now()} -- Running contrast check...')
            contrast_check = [True if x <= 5 else False for x in contrast]
            print(f'Current time: {datetime.datetime.now()} -- Finished contrast check')

            # Split the parameters into chunks to pass to workers
            print(f'Current time: {datetime.datetime.now()} -- Splitting parameters into chunks to pass to workers...')
            divisor = num_generated // cpu_ct
            n_divisor = num_generated // divisor

            # Create an array with all necessary arguments for processing
            param_stack = [[period[j],ecc[j],
                           arg_peri[j],phase[j],self.MJD]
                           for j in range(num_generated)]
            # Move into parallel processing
            #  Create Processes
            pool = mp.Pool(cpu_ct)

            # Use Pool to calculate RVs
            print(f'Current time: {datetime.datetime.now()} -- Parallel processes are now calculating RVs (passed to workers!)...')

            result = pool.starmap(self.calculate_unscaled_RV, param_stack, chunksize=divisor)
            print(f'Current time: {datetime.datetime.now()} -- Parallel processes done calculating RVs!')
            unscaled_rv = np.asarray(result)

            # The barycentric semimajor axes are related by the mass ratio
            a_prim = a * (mass_ratio / (mass_ratio + 1))
            a_cmp = a - a_prim
            
            # Most of the terms in K are constants for a given orbit, 
            # then rescale by individual a values
            sin_i = np.sin(np.arccos(cos_i))
            K_term = (2*np.pi/period) * sin_i / np.sqrt(1-ecc**2)
            prim_K = K_term * a_prim
            cmp_K = K_term * a_cmp

            # scaled rvs equal unscaled rv * K
            # and the companion moves opposite to the primary
            urv_t = unscaled_rv.T
            prim_rv = (prim_K * urv_t).T
            cmp_rv = -1 * (cmp_K * urv_t).T

            max_delta_rv = np.max(np.absolute(np.subtract(prim_rv, cmp_rv)), axis=1)
            print(f'Current time: {datetime.datetime.now()} -- Finished concatenating RV results...')


            # Determine the overall predicted RV
#            comb_rvs = np.hstack([prim_rv,cmp_rv])
#            prim_lum_arr = np.full_like(prim_rv,fill_value=prim_lum)
#            cmp_lum_arr = np.repeat(cmp_lum,nobs).reshape((-1,nobs))
#            comb_weights = np.hstack([prim_lum_arr, cmp_lum_arr])
#            self.predicted_RV = np.average(comb_rvs,axis=1,weights=comb_weights)
            self.predicted_RV = np.zeros_like(prim_rv)

            # SB1s have high contrast
            sb1 = contrast > 5
            self.predicted_RV[sb1,:] = prim_rv[sb1,:]

            sb2 = np.where(contrast>=5)[0]
            for i in sb2:
                self.predicted_RV[i] = np.average([prim_rv[i],cmp_rv[i]],
                                                axis=0,
                                                weights=[prim_lum,cmp_lum[i]])



            # Use Pool to calculate zero point
            print(f'Current time: {datetime.datetime.now()} -- Calculating zero point...')
            zp_params = [[self.experimental_RV, self.measurement_error, self.predicted_RV[j]] for j in range(num_generated)]

            zero_points = pool.starmap(zero_point_fit_parallel, zp_params,
                                       chunksize=divisor)
            print(f'Current time: {datetime.datetime.now()} -- Calculated zero point!')

            # Shift all by zero point
            print(f'Current time: {datetime.datetime.now()} -- Shifting all by zero point...')
            self.predicted_RV = self.predicted_RV + zero_points
            print(f'Current time: {datetime.datetime.now()} -- Shifted all by zero point!')

            # Compare experimental and predicted RVs
            print(f'Current time: {datetime.datetime.now()} -- Comparing experimental and predicted RVs...')
            
            # Calculate amp
            amp = np.ptp(self.predicted_RV, axis=1)
            

            # Calculate chi^2            
            chi_sq_numer = np.square(np.subtract(self.experimental_RV, self.predicted_RV))
            chi_sq_denom = self.measurement_error**2 + self.added_jitter**2
            chi_squared0 = chi_sq_numer / chi_sq_denom
            chi_squared = np.sum(chi_squared0,axis=1)

            # Calculate prob
            ndatesm1 = len(self.MJD)-1
            # TODO: can this be further streamlined?
            pr_params = [[chi2,ndatesm1] for chi2 in chi_squared]
            # prob = [stats.chi2.cdf(chi_squared[i], ndatesm1) for i in range(num_generated)]
            pr_results = pool.starmap(stats.chi2.cdf,pr_params,
                                      chunksize=divisor)
            prob = np.asarray(pr_results)
            pool.close()


        # End Parallelized

        # Reject things with a rejection probability greater than 0.997, corresponding to 3 sigma
        rv_fit_reject = np.random.rand(num_generated) < prob

        # Check amplitude and resolution
        above_amplitude = np.abs(amp) > self.rv_floor
        self.amp = amp
        visible_sb2 = (contrast <= 5) & (max_delta_rv > delta_v)
        invisible_sb2 = (contrast <= 5) & (max_delta_rv < delta_v)
        self.b_type = np.full(num_generated,"","U14")
        self.b_type[contrast > 5] = "SB1"
        self.b_type[visible_sb2] = "Resolved SB2"
        self.b_type[invisible_sb2] = "Unresolved SB2"

        self.rv_reject_list = (rv_fit_reject & above_amplitude) | visible_sb2
        
        return self.rv_reject_list

    def calculate_jitter(self, predicted, experimental, error, threshold=0.00001):
        # Calculates the stellar jitter term necessary to produce a chi_squared value of 1
        # Inputs: Measured RV, Predicted RV, Measurement Error
        # Outputs: Jitter
        print(f'Current time: {datetime.datetime.now()} -- Calculating jitter')
        a = [experimental[i]-predicted[i] for i in range(0, len(experimental))]
        b = np.std(a)**2
        c = np.median(error)**2
        initial_guess = np.sqrt(b+c)
        prev_jitter = 0
        new_jitter = initial_guess
        itr = 0
        # Using Newton's method to iterate
        while abs(new_jitter - prev_jitter) > threshold and itr < 10000:
            # Change prev_jitter
            prev_jitter = new_jitter
            # Calculate f(j) and f'(j)
            chi_square = sum([(experimental[i] - predicted[i]) ** 2 / (error[i]**2 + prev_jitter**2) for i in range(0, len(error))])
            chi_square_prime = sum([-(experimental[i] - predicted[i])**2 / (error[i]**2 + prev_jitter**2)**2 * 2 * prev_jitter for i in range(0, len(error))])
            # Calculate new jitter
            new_jitter = prev_jitter - (chi_square - 1) / chi_square_prime
            itr += 1
        if abs(new_jitter - prev_jitter) < threshold:
            print(f'Current time: {datetime.datetime.now()} -- Finished calculating jitter')
            return new_jitter
        else:
            print(f'Current time: {datetime.datetime.now()} -- Finished calculating jitter')
            return -1

    @ staticmethod
    def calculate_RV(period, mass_ratio, a, e, cos_i, arg_peri, phase, MJD):
        # Calculates the RVs for each item when passed arrays of orbital parameters
        # Inputs: Arrays of Period, Mass Ratio, Semi-Major Axis, eccentricity, inclination, arg peri, phase, calculation times
        # Outputs: Velocity Semi-Amplitude (km/s), RVs at each time in MJD
        print(f'\n\nCurrent time: {datetime.datetime.now()} -- Calculating RVs {mp.current_process()}...')
        
        sin_i = np.sin(np.arccos(cos_i))


        n = len(period)
        ndates = len(MJD)
        # Create a blank array with one row for every input period
        # and each row contains RVs equivalent to the observation dates
        RV = np.zeros((n,ndates))
        # RV = [[0.0 for i in range(len(MJD))] for j in range(n)]
        
        a_star = a * (mass_ratio / (mass_ratio + 1))
        # a_star = np.multiply(a, np.divide(mass_ratio, np.add(mass_ratio, 1)))

        K = (2*np.pi/period) * (a_star * sin_i) / np.sqrt(1-ecc**2)
        # K = np.multiply(np.divide((2 * np.pi), period),np.divide(np.multiply(a_star, sin_i), np.sqrt((1 - np.square(ecc)))))  # AU/days
        K = K * 1731.48  # km/s

        print(f'Current time: {datetime.datetime.now()} -- Iterating over companions {mp.current_process()}...')
        for i in range(n):  # Iterate over companions
            for j in range(0, len(MJD)):  # Iterate over times
                # Find E
                M = 2 * np.pi * MJD[j] / period[i] - phase[i]
                prev_E = 0.0
                current_E = M

                while abs(current_E - prev_E) > 0.00001:
                    prev_E = current_E
                    current_E = M + e[i] * np.sin(prev_E)

                # Find true anomaly, f
                f = 2 * np.arctan2(np.tan(current_E / 2), np.sqrt((1 - e[i]) / (1 + e[i])))
                # Find predicted RV
                RV[i][j] = K[i] * (np.sin(arg_peri[i] + f) + e[i] * np.sin(arg_peri[i]))  # km/s
        print(f'Current time: {datetime.datetime.now()} -- Finished iterating over companions {mp.current_process()}...')


        print(f'Current time: {datetime.datetime.now()} -- Finished calculating RVs! {mp.current_process()}')
        return K, RV

    @ staticmethod
    def calculate_unscaled_RV(period, ecc, arg_peri, phase, MJD):
        """
        Calculate the *unscaled* radial velocity for a given period,
        eccentricity, arg_periastron, and phase at the given MJD dates.
        Arrays are NOT accepted for any argument except MJD. 

        Returns only the unscaled RV. 
        """
        
        ndates = len(MJD)
        # RV = np.zeros(ndates)
        
        f_ecc_term = np.sqrt((1 - ecc) / (1 + ecc))
        pi2per = 2*np.pi/period

        true_anomaly = np.zeros(ndates)

        for j in range(0, len(MJD)):  # Iterate over times
            # Find E
            M = pi2per * MJD[j] - phase
            prev_E = 0.0
            current_E = M

            while abs(current_E - prev_E) > 0.00001:
                prev_E = current_E
                current_E = M + ecc * np.sin(prev_E)

            # Find true anomaly, f
            true_anomaly[j] = 2 * np.arctan2(np.tan(current_E / 2), f_ecc_term)
            # Find predicted RV
        RV = np.sin(arg_peri + true_anomaly) + ecc * np.sin(arg_peri)  # km/s

        return RV

    def read_in_rv(self):

        if os.path.exists(self.rv_filename)==False:
            return -31

        try:
            rv = Table.read(self.rv_filename, format='ascii', delimiter=' ', fast_reader=False)
            col_names = list(rv.columns)
            assert len(col_names) == 3
        except:
            try:
                rv = Table.read(self.rv_filename, format='ascii', delimiter='\t', fast_reader=False)
            except Exception as e:
                return -32
        
        print(f'Current time: {datetime.datetime.now()} -- Renaming and sorting RV columns...')
        # rename the columns so they match what i want
        col_names = list(rv.columns)
        if len(col_names) != 3:
            return -32
        rv.rename_column(col_names[0],'JD')
        rv.rename_column(col_names[1],'RV')
        rv.rename_column(col_names[2],'RVerr')

        rv.sort('JD')

        t_0 = 2400000.5
        self.MJD = rv['JD']-t_0
        self.experimental_RV = rv['RV']
        self.measurement_error = rv['RVerr']
        print(f'Current time: {datetime.datetime.now()} -- Finished renaming and sorting RV columns')

        return 0

    def restore_defaults(self):
        self.rv_filename = ''
        self.added_jitter = 0.02  # km/s
        self.rv_floor = 0.02  # km/s
        self.rv_reject_list = []
        self.jitter_reject_list = []
        self.MJD = []
        self.experimental_RV = []
        self.measurement_error = []
        self.predicted_RV = []
        self.jitter = []
        gc.collect()

    def zero_point_model(self, prediction, a):
        # Alters the prediction by some zero point shift, a
        y = prediction + a
        return y

    def load_stellar_model(self, filter, star_age):
        print(f'Current time: {datetime.datetime.now()} -- Loading stellar model from the actual function itself...')
        # Read in file containing stellar model with the filter needed
        # RV always loads G filter
        print(f'Current time: {datetime.datetime.now()} -- Reading file containing stellar model with necessary filter...')
        if filter == 'J' or filter == 'H' or filter == 'K':  # 2MASS filters
            print(f'Current time: {datetime.datetime.now()} -- Reading BHAC file and creating table...')
            model_chart = {}
            BHAC_file = f'{os.path.join(repo_path, "reference_data/BHAC15_2MASS.txt")}'
            with open(BHAC_file, 'r') as content_file:
                content = content_file.read()
            tables = content.split(
                sep='\n-----------------------------------------------------------------------------------------------\n')
            tables = [x for x in tables if len(x) > 1]
            print("\nTables: {tables}")
            print(f'Current time: {datetime.datetime.now()} -- Finished reading BHAC file and creating table')

            print(f'Current time: {datetime.datetime.now()} -- ...')
            for table in tables:
                n = table.find('M')
                time_segment = table[0:n]
                # print(f"\nTime segment: {time_segment}")
                table_segment = table[n:]
                # print(f"\nTable segment: {time_segment}")
                age = float(time_segment[time_segment.find('=') + 1:])
                year_chart = Table.read(table_segment, format='ascii', fast_reader=False)
                if filter == 'J':
                    year_chart = year_chart['M/Ms', 'Mj']
                    year_chart.rename_column('Mj', 'Mag')
                elif filter == 'H':
                    year_chart = year_chart['M/Ms', 'Mh']
                    year_chart.rename_column('Mh', 'Mag')
                elif filter == 'K':
                    year_chart = year_chart['M/Ms', 'Mk']
                    year_chart.rename_column('Mk', 'Mag')
                # year_chart.show_in_browser(jsviewer=True, tableid="year_chart")
                model_chart[age] = year_chart
                # print(f"\nModel chart[age]: {model_chart[age]}")

        elif filter == 'G' or filter == 'R' or filter == 'I':
            model_chart = {}
            BHAC_file = f'{os.path.join(repo_path, "reference_data/BHAC15_CFHT.txt")}'
            print(f'Current time: {datetime.datetime.now()} -- Reading BHAC file and creating table...')
            with open(BHAC_file, 'r') as content_file:
                content = content_file.read()
            tables = content.split(
                sep='\n----------------------------------------------------------------------------------------------------------------------------------------\n')
            tables = [x for x in tables if len(x) > 1]
            print(f'Current time: {datetime.datetime.now()} -- Finished reading BHAC file and creating table...')

            for table in tables:
                n = table.find('M')
                time_segment = table[0:n]
                table_segment = table[n:]

                age = float(time_segment[time_segment.find('=') + 1:])
                year_chart = Table.read(table_segment, format='ascii', fast_reader=False)
                if filter == 'G':
                    year_chart = year_chart['M/Ms', 'G', 'L/Ls']
                    year_chart.rename_column('G', 'Mag')
                elif filter == 'R':
                    year_chart = year_chart['M/Ms', 'R', 'L/Ls']
                    year_chart.rename_column('R', 'Mag')
                elif filter == 'I':
                    year_chart = year_chart['M/Ms', 'I', 'L/Ls']
                    year_chart.rename_column('I', 'Mag')
                model_chart[age] = year_chart
           
            # print(f"\nTime segment: {time_segment}")
            # print(f"\nTable segment: {time_segment}")
            # year_chart.show_in_browser(jsviewer=True, tableid="year_chart")
            # print(f"\nModel chart[age]: {model_chart[age]}")

        ages = np.array(list(model_chart.keys()))
        #  Check if age is modeled, and if it is simply return that table
        if star_age in ages:
            return model_chart[star_age]
        # If the age is not included in the models, linearly interpolate the parameters from included ages
        #  Find ages above and below the desired age
        print(f'Current time: {datetime.datetime.now()} -- Find ages above and below desired age...')
        diff = [star_age - x for x in ages]
        low_age = ages[np.argmin([x if x > 0 else float('inf') for x in diff])]
        high_age = ages[np.argmin([abs(x) if x < 0 else float('inf') for x in diff])]
        young_chart = model_chart[low_age]
        old_chart = model_chart[high_age]
        print(f'Current time: {datetime.datetime.now()} -- Finished find ages above and below desired age')


        #  Get masses
        print(f'Current time: {datetime.datetime.now()} -- Getting masses...')
        common_mass = np.intersect1d(young_chart['M/Ms'], old_chart['M/Ms'])

        young_chart = young_chart[[x in common_mass for x in  young_chart['M/Ms']]]
        old_chart = old_chart[[x in common_mass for x in old_chart['M/Ms']]]

        new_model = Table()
        new_model['M/Ms'] = common_mass

        #  Interpolate
        print(f'Current time: {datetime.datetime.now()} -- Interpolating masses...')
        for col in model_chart[low_age].colnames[1:]:
            col_list = []
            for i in range(len(common_mass)):
                # interpolate the new value of the parameter for each mass
                xs = [low_age, high_age]
                ys = [young_chart[col][i], old_chart[col][i]]
                f = scipy.interpolate.interp1d(xs, ys, kind='linear')
                col_list.append(f(star_age))
            # add to table
            new_model[col] = col_list
            print(f'Current time: {datetime.datetime.now()} -- Finished interpolating masses...')
        print(f'Current time: {datetime.datetime.now()} -- Finished getting masses...')

        return new_model
    
def apply_ptp(a=[0,1,2]):
    return np.ptp(a)
