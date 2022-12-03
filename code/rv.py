from datetime import datetime as dt
import datetime
import numpy as np
import scipy as scipy
import scipy.stats as stats
from astropy.table import Table
from time import time
import multiprocessing as mp
import gc
import warnings
import os
from astropy.utils.exceptions import AstropyWarning
import logging
from timeit import timeit
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

today = dt.today().isoformat().split("T")[0]
global repo_path
repo_path = os.getenv('MOLOC').replace("\\", "/")

def calculate_RV_parallel(period, mass_ratio, a, e, cos_i, arg_peri, phase, MJD, calc):
    # Exactly the same as calculate_RV, but with an extra parameter stating whether you need to calculate RV
    # Calculates the RVs for each item when passed arrays of orbital parameters
    # Inputs: Arrays of Period, Mass Ratio, Semi-Major Axis, eccentricity, inclination, arg peri, phase, calculation times
    # Outputs: Velocity Semi-Amplitude (km/s), RVs at each time in MJD


    # print(f'Current time: {datetime.datetime.now()} -- Calculating RV Para
    sin_i = np.sin(np.arccos(cos_i))

    n = len(period)
    RV = [[0.0 for i in range(len(MJD))] for j in range(n)]
    a_star = np.multiply(a, np.divide(mass_ratio, np.add(mass_ratio, 1)))
    K = np.multiply(np.divide((2 * np.pi), period),np.divide(np.multiply(a_star, sin_i), np.sqrt((1 - np.square(e)))))  # AU/days
    K = np.multiply(K, 1731.48)  # km/s

    for i in range(n):  # Iterate over companions
        if calc[i]:
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

    return K, RV

def zero_point_model(prediction, a):
    # Alters the prediction by some zero point shift, a
    y = prediction + a
    return y

def zero_point_fit_parallel(experimental, errors, predicted):

    A = [0.] * len(predicted)
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
        e = self.companions.ecc
        cos_i = self.companions.cos_i
        arg_peri = self.companions.arg_peri
        phase = self.companions.phase
        print(f'Current time: {datetime.datetime.now()} -- Finished unpacking parameters')


        # Determine velocity limit
        print(f'Current time: {datetime.datetime.now()} -- Determining velocity limit...')
        delta_v = 2.998e5 / self.resolution  # m/s

        t = time()

        # Determine contrast
        cmp_mass = np.multiply(self.star_mass, mass_ratio)  # companion mass in solar masses

        f_mag = scipy.interpolate.interp1d(self.model['M/Ms'], self.model['Mag'], kind='cubic', fill_value='extrapolate')
        # self.model.show_in_browser(jsviewer=True, tableid="self.model")


        prim_model_mag = f_mag(self.star_mass)
        cmp_model_mag = [f_mag(x) if x >= self.model['M/Ms'][0] else float('inf') for x in cmp_mass]  # companion mags, assign infinite magnitude if below lowest modeled mass
        contrast = np.subtract(cmp_model_mag, prim_model_mag)
        print(f'Current time: {datetime.datetime.now()} -- Finished determining velocity limit')


        # Determine luminosity
        print(f'Current time: {datetime.datetime.now()} -- Determining luminosity...')
        f_lum = scipy.interpolate.interp1d(self.model['M/Ms'], self.model['L/Ls'], kind='cubic', fill_value='extrapolate')
        prim_lum = np.power(10, f_lum(self.star_mass))  # primary luminosity
        cmp_lum = [np.power(10, f_lum(x)) if x >= self.model['M/Ms'][0] else 0. for x in cmp_mass]  # companion luminosity, 0 if below limit
        print(f'Current time: {datetime.datetime.now()} -- Determined luminosity')

        # Choose to run either parallelized or non-parallelized based on the number of generated companions
        print(f'Current time: {datetime.datetime.now()} -- Determining whether or not to run parallelized...')
        if num_generated < 50000:  # Serial
            print(f'Current time: {datetime.datetime.now()} -- Decided to run non-parallelzied')
            self.predicted_RV = [np.zeros(len(self.MJD)) for x in range(num_generated)] # pre-allocate

            # calculate RV curves
            prim_K, prim_rv = self.calculate_RV(period, mass_ratio, a, e, cos_i, arg_peri, phase, self.MJD)
            cmp_K, cmp_rv = self.calculate_RV(period, np.divide(1, mass_ratio), a, e, cos_i, arg_peri, phase, self.MJD)
            cmp_rv = np.multiply(-1, cmp_rv)

            max_delta_rv = np.max(np.absolute(np.subtract(prim_rv, cmp_rv)), axis=1)

            # Determine the overall predicted RV
            print(f'Current time: {datetime.datetime.now()} -- Determining overall predicted RV...')
            for i in range(num_generated):
                if contrast[i] > 5:
                    # SB1
                    self.predicted_RV[i] = prim_rv[i]
                elif contrast[i] <= 5:
                    # SB2, looks like SB1. RV is weighted average
                    rv = np.average([prim_rv[i], cmp_rv[i]], axis=0, weights=[prim_lum, cmp_lum[i]])
                    self.predicted_RV[i] = rv

            for i in range(0, num_generated):
                # Fit the zero point
                [zero_point], pcov = scipy.optimize.curve_fit(self.zero_point_model, self.predicted_RV[i], self.experimental_RV, sigma=self.measurement_error)
                # Shift all predicted values by the zero_point
                self.predicted_RV[i] += zero_point
            print(f'Current time: {datetime.datetime.now()} -- Determined overall predicted RV')

            print(f'Current time: {datetime.datetime.now()} -- Comparing experimental RV to predicted RV...')
            # Compare experimental and predicted RVs
            amp = [np.ptp(self.predicted_RV[i]) for i in range(num_generated)]
            chi_squared = [sum(np.divide(np.square(np.subtract(self.experimental_RV, self.predicted_RV[i])),
                        np.add(np.square(self.measurement_error), self.added_jitter ** 2))) for i in range(num_generated)]
            # The degrees of freedom is equal to (N-1)+1, for the number of data points and the applied velocity shift
            prob = [stats.chi2.cdf(chi_squared[i], len(self.MJD)-1) for i in range(0, num_generated)]
            print(f'Current time: {datetime.datetime.now()} -- Finished comparing experimental RV to predicted RV')

        else:  
            # Parallel
            print(f'Current time: {datetime.datetime.now()} -- Decided to run parallelzied')
            # Determine cpu count
            try:
                cpu_ct = len(os.sched_getaffinity(0))-1
            except AttributeError:
                cpu_ct = mp.cpu_count()-1

            print(f'Current time: {datetime.datetime.now()} -- Running contrast check...')
            contrast_check = [True if x <= 5 else False for x in contrast]
            print(f'Current time: {datetime.datetime.now()} -- Finished contrast check')

            # Split the parameters into chunks to pass to workers
            print(f'Current time: {datetime.datetime.now()} -- Splitting parameters into chunks to pass to workers...')
            divisor = int(np.ceil(min(num_generated / cpu_ct, 200000)))
            n_divisor = int(np.ceil(num_generated / divisor))

            # Parameters are split into chunks of size divisor
            period = [period[i:i + divisor] for i in range(0, num_generated, divisor)]
            mass_ratio = [mass_ratio[i:i + divisor] for i in range(0, num_generated, divisor)]
            a = [a[i:i + divisor] for i in range(0, num_generated, divisor)]
            e = [e[i:i + divisor] for i in range(0, num_generated, divisor)]
            cos_i = [cos_i[i:i + divisor] for i in range(0, num_generated, divisor)]
            arg_peri = [arg_peri[i:i + divisor] for i in range(0, num_generated, divisor)]
            phase = [phase[i:i + divisor] for i in range(0, num_generated, divisor)]
            contrast_check = [contrast_check[i:i+divisor] for i  in range(0,  num_generated, divisor)]
            print(f'Current time: {datetime.datetime.now()} -- Parameters have been split...')

            # Move into parallel processing
            #  Create Processes
            print(f'Current time: {datetime.datetime.now()} -- Creating parallel processes...')
            pool = mp.Pool(cpu_ct)
            print(f'Current time: {datetime.datetime.now()} -- Parallel processes have been created')

            # Use Pool to calculate RVs
            print(f'Current time: {datetime.datetime.now()} -- Parallel processes are now calculating RVs (passed to workers!)...')
            prim_results = pool.starmap(self.calculate_RV, [(period[j], mass_ratio[j], a[j], e[j], cos_i[j], arg_peri[j], phase[j], self.MJD) for j in range(n_divisor)])
            cmp_results = pool.starmap(calculate_RV_parallel, [(period[j], np.divide(1, mass_ratio[j]), a[j],
                    e[j], cos_i[j], arg_peri[j], phase[j], self.MJD, contrast_check[j]) for j in range(n_divisor)])
            print(f'Current time: {datetime.datetime.now()} -- Parallel processes done calculating RVs!')

            # Concatenate Results
            print(f'Current time: {datetime.datetime.now()} -- Concatenating RV results...')
            prim_K = np.hstack([prim_results[i][0] for i in range(int(np.ceil(num_generated / divisor)))])
            prim_rv = np.vstack([prim_results[i][1] for i in range(int(np.ceil(num_generated / divisor)))])
            cmp_K = np.hstack([cmp_results[i][0] for i in range(int(np.ceil(num_generated / divisor)))])
            cmp_rv = np.vstack([cmp_results[i][1] for i in range(int(np.ceil(num_generated / divisor)))])
            cmp_rv = np.multiply(-1, cmp_rv)

            max_delta_rv = np.max(np.absolute(np.subtract(prim_rv, cmp_rv)), axis=1)
            print(f'Current time: {datetime.datetime.now()} -- Finished concatenating RV results...')


            # Determine the overall predicted RV
            print(f'Current time: {datetime.datetime.now()} -- Determining overall predicted RV...')
            self.predicted_RV = [np.zeros(len(self.MJD)) for x in range(num_generated)]  # pre-allocate

            for i in range(num_generated):
                if contrast[i] > 5:
                    # SB1
                    self.predicted_RV[i] = prim_rv[i]
                elif contrast[i] <= 5:
                    # SB2, looks like SB1. RV is weighted average
                    rv = np.average([prim_rv[i], cmp_rv[i]], axis=0, weights=[prim_lum, cmp_lum[i]])
                    self.predicted_RV[i] = rv
                    
            print(f'Current time: {datetime.datetime.now()} -- Determined overall predicted RV')

            # Use Pool to calculate zero point
            print(f'Current time: {datetime.datetime.now()} -- Calculating zero point...')
            split_RV  = [self.predicted_RV[i:i+divisor] for i  in range(0,  num_generated, divisor)]
            zero_points = pool.starmap(zero_point_fit_parallel, [(self.experimental_RV, self.measurement_error, split_RV[j]) for j in range(n_divisor)])
            zero_points = np.concatenate(zero_points, axis=0)
            pool.close()
            print(f'Current time: {datetime.datetime.now()} -- Calculated zero point!')

            # Shift all by zero point
            print(f'Current time: {datetime.datetime.now()} -- Shifting all by zero point...')
            self.predicted_RV = [np.add(self.predicted_RV[i], zero_points[i]) for i in range(num_generated)]
            print(f'Current time: {datetime.datetime.now()} -- Shifted all by zero point!')

            # Compare experimental and predicted RVs
            print(f'Current time: {datetime.datetime.now()} -- Comparing experimental and predicted RVs...')
            # self.predicted_RV.show_in_browser(js_viewer=True)
            # print(f"---------------------------------------------Predicted RV: {type(self.predicted_RV)}")
            # print(f"---------------------------------------------Predicted RV length: {len(self.predicted_RV)}")
            # print(f"---------------------------------------------Num generated: {num_generated}")
            # print(f"---------------------------------------------Predicted RV[num]: {self.predicted_RV[range(num_generated-1)]}")
            # print(f"---------------------------------------------Predicted RV[i]:")
            # for i in range(3):
            #     print(self.predicted_RV[i])
            
            print(f'Current time: {datetime.datetime.now()} ----------------------------------- Pre map')
            amp_test = list(map(apply_ptp, self.predicted_RV))
            print(f'Current time: {datetime.datetime.now()} ----------------------------------- Post map')

            print(f'Current time: {datetime.datetime.now()} ----------------------------------- Pre amp loop')
            amp = [np.ptp(self.predicted_RV[i]) for i in range(num_generated)]
            print(f'Current time: {datetime.datetime.now()} ----------------------------------- Post amp loop')

            # if amp_test == amp:
            #     print("They are equal!")
            #     print(len(amp_test))
            #     print(len(amp))
            # else:
            #     print("Not equal :(((")
            #     print(len(amp_test))
            #     print(len(amp))

            # amp = [np.ptp(self.predicted_RV[range(num_generated-1)])]
            chi_squared = [sum(np.divide(np.square(np.subtract(self.experimental_RV, self.predicted_RV[i])),
                           np.add(np.square(self.measurement_error), self.added_jitter**2))) for i in range(num_generated)]
            prob = [stats.chi2.cdf(chi_squared[i], len(self.MJD)-1) for i in range(0, num_generated)]
            
            print(f'Current time: {datetime.datetime.now()} -- Compared experimental and predicted RVs!')
        # End Parallelized

        # Reject things with a rejection probability greater than 0.997, corresponding to 3 sigma
        print(f'Current time: {datetime.datetime.now()} -- Rejecting unlikely companions...')
        rv_fit_reject = np.array([True if np.random.rand() < x else False for x in prob])
        # Check amplitude and resolution
        print(f'Current time: {datetime.datetime.now()} -- Checking amplitude and resolution...')
        above_amplitude = np.array([True if abs(x) > self.rv_floor else False for x in amp])
        self.amp = amp
        visible_sb2 = np.array([True if contrast[i] < 5 and max_delta_rv[i] > delta_v else False for i in range(num_generated)])
        self.b_type = ['              ']*num_generated
        for i in range(num_generated):
            if contrast[i] > 5:
                self.b_type[i] = 'SB1'
            elif contrast[i] <=5 and max_delta_rv[i] < delta_v:
                self.b_type[i] = 'Unresolved SB2'
            elif contrast[i] <=5 and max_delta_rv[i] > delta_v:
                self.b_type[i] = 'Resolved SB2'
        print(f'Current time: {datetime.datetime.now()} -- Finished checking amplitude and resolution')

        self.rv_reject_list = np.array([True if (rv_fit_reject[i] and above_amplitude[i]) or visible_sb2[i] else False for i in range(num_generated)])
        print(f'Current time: {datetime.datetime.now()} -- Rejected unlikely companions')
        
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
        print(f'Current time: {datetime.datetime.now()} -- Calculating RVs {mp.current_process()}...')

        sin_i = np.sin(np.arccos(cos_i))

        n = len(period)
        RV = [[0.0 for i in range(len(MJD))] for j in range(n)]
        a_star = np.multiply(a, np.divide(mass_ratio, np.add(mass_ratio, 1)))
        K = np.multiply(np.divide((2 * np.pi), period),np.divide(np.multiply(a_star, sin_i), np.sqrt((1 - np.square(e)))))  # AU/days
        K = np.multiply(K, 1731.48)  # km/s
        
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
      
        print(f'Current time: {datetime.datetime.now()} -- Finished calculating RVs! {mp.current_process()}')
        return K, RV

    def read_in_rv(self):
        print(f'Current time: {datetime.datetime.now()} -- Reading in RVs...')

        try:
            rv_in = open(self.rv_filename, 'r')
            print(f'Current time: {datetime.datetime.now()} -- Finished reading in RVs!')
        except FileNotFoundError:
            return -31

        try:
            rv = Table.read(self.rv_filename, format='ascii', delimiter=' ', fast_reader=False)
            col_names = list(rv.columns)
            assert len(col_names) == 3
            print(f'Current time: {datetime.datetime.now()} -- Finished reading in RVs!')
        except:
            try:
                rv = Table.read(self.rv_filename, format='ascii', delimiter='\t', fast_reader=False)
                print(f'Current time: {datetime.datetime.now()} -- Finished reading in RVs!')
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
        print(f'Current time: {datetime.datetime.now()} -- Running zero point model...')
        y = prediction + a
        print(f'Current time: {datetime.datetime.now()} -- Finished running zero point model')
        return y

    def load_stellar_model(self, filter, star_age):
        print(f'Current time: {datetime.datetime.now()} -- Loading stellar model from the actual function itself...')
        # Read in file containing stellar model with the filter needed
        # RV always loads G filter
        print(f'Current time: {datetime.datetime.now()} -- Reading file containing stellar model with necessary filter...')
        if filter == 'J' or filter == 'H' or filter == 'K':  # 2MASS filters
            print(f'Current time: {datetime.datetime.now()} -- Reading BHAC file and creating table...')
            model_chart = {}
            BHAC_file = f'{os.path.join(repo_path, "code/BHAC15_2MASS.txt")}'
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
            BHAC_file = f'{os.path.join(repo_path, "code/BHAC15_CFHT.txt")}'
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