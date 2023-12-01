from datetime import datetime
import warnings
import logging

import h5py
import numpy as np
import scipy as scipy
import scipy.stats as stats
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

today = datetime.today().isoformat().split("T")[0]

class Companions:
    # class variables
    __P = []  # days
    __mass_ratio = []  # fraction
    __cos_i = []  # unitless
    __a = []  # AU
    __ecc = []  # unitless
    __arg_peri = []  # radians
    __phase = []  # radians
    num_generated = 0

    # class functions
    def __init__(self, num_generated, limits, star_mass, pd_mu, pd_sig, q_exp):
        self.num_generated = num_generated
        self.limits = limits
        self.star_mass = star_mass
        self.mu_log_P = pd_mu
        self.sig_log_P = pd_sig
        self.q_exp = q_exp

    def generate(self):

        np.random.seed(999)

        # If the user has fixed P, a, AND mass_ratio, then we can't continue
        if ((self.limits[0] is not None) and (self.limits[15] is not None) and
            (self.limits[12] is not None)):
            raise ValueError("You can't fix P, a, and q at the same time")

        self.generate_paq()
        self.generate_inclination()
        self.generate_phase()
        self.generate_eccentricity()
        self.generate_arg_peri()

    def generate_paq(self):
        """
        generate Period (days), Mass Ratio and Semi-Major Axis (AU)

        """
        G = 39.478  # Gravitational constant in AU^3/years^2*M_solar
        P_fixed = self.limits[0]
        P_lower = self.limits[1]
        P_upper = self.limits[2]
        mass_fixed = self.limits[12]
        mass_lower = self.limits[13]
        mass_upper = self.limits[14]
        a_fixed = self.limits[15]
        a_lower = self.limits[16]
        a_upper = self.limits[17]

        print('Generation Mass Ratio Exponent: ', self.q_exp)
        dq = 1e-4  # coarseness

        # Default lower limit on P
        if P_lower is None or P_lower < 0.1:
            P_lower = 0.1

        if mass_fixed is not None:
            self.mass_ratio = np.full(self.num_generated,fill_value=mass_fixed)

        if a_fixed is not None:
            self.a = np.full(self.num_generated,fill_value=a_fixed)

        # We will always generate P. Either it's fixed, OR the
        # other two values are fixed, OR we generate it independently
        if P_fixed is not None:
            self.P = np.full(self.num_generated,fill_value=P_fixed)
        # with a and q, calculate P directly
        elif a_fixed is not None and mass_fixed is not None:
            self.P = np.sqrt( (4 * np.pi ** 2 * self.a ** 3) / (G * self.star_mass * (1 + self.mass_ratio))) * 365
        else:
            tn_low = (P_lower-self.mu_log_P)/self.sig_log_P
            if P_upper is None:
                tn_up = np.inf
            else:
                tn_up = (P_upper-self.mu_log_P)/self.sig_log_P
            P_tn = stats.truncnorm(tn_low,tn_up,
                                   loc=self.mu_log_P,scale=self.sig_log_P)
            log_P = P_tn.rvs(self.num_generated)
            self.P = 10**log_P

        # Now we handle all other cases, either with two fixed parameters
        # including P, or with only one fixed parameter (or none)

        # With a and P, calculate q directly
        # If a is fixed, then P and q cannot be generated independently
        # Either way, we already generate P above
        # and now we will calculate the mass ratio accordingly
        if a_fixed is not None and mass_fixed is None:
            self.mass_ratio = np.array(4*np.pi**2 * self.a**3 / (G * self.star_mass * (self.P / 365)** 2) - 1)
            if np.any(self.mass_ratio > 1):
                # The fixed values of period and semi-major axis require a mass ratio greater than one
                return -12

        # We've handled all cases where values can be fixed. 
        # Now we generate the mass ratio randomly if needed, 
        # and calculate a if it wasn't fixed
        else:
            # We already handled this case above
            if mass_fixed is not None:
                pass

            # Flat IMF
            elif self.q_exp == 0.0:
                # This generates all the  way to q=0 for a uniform distribution
                if mass_lower is None:
                    mass_lower = 0.
                if mass_upper is None:
                    mass_upper = 1.
                self.mass_ratio = np.random.uniform(mass_lower, mass_upper, self.num_generated)

            # IMF with an input slope
            else:
                # This applies a lower mass limit for the power law distribution
                if mass_lower is None:
                    mass_lower = 0.01
                if mass_upper is None:
                    mass_upper = 1.
                q_range = np.arange(mass_lower, mass_upper + dq, dq)  # range of possible q
                p = np.power(q_range, self.q_exp)    # probabilities
                p = p/np.sum(p) # normalized probabilities
                self.mass_ratio = np.random.choice(q_range, p=p, size=self.num_generated)
            # Now calculate Semi-Major Axis using P and m
            self.a = ((self.P / 365)**2 * G * self.star_mass*(1 + self.mass_ratio)/(4 * np.pi ** 2))**(1/3)  # AU         
            

    def generate_inclination(self):
        # Inclination
        cos_i_fixed = self.limits[3]
        cos_i_lower = self.limits[4]
        cos_i_upper = self.limits[5]

        # todo figure out what to do about stellar radius in this calculation
        if cos_i_lower == 'transit' or cos_i_upper == 'transit' or cos_i_fixed == 'transit':
            self.cos_i = np.zeros(self.num_generated)
            b_limit = 1.01
            R_star = 1.  # R_sun
            r_sun_to_au = 0.00465
            for i in range(0, self.num_generated):
                # Determine the upper limit on cos(i), based on semi-major axis, stellar radius, and impact parameter
                cos_i_upper = np.min([(b_limit * r_sun_to_au * R_star)/(self.a[i]), 1.0])
                #  Generate
                self.cos_i[i] = np.random.uniform(0., cos_i_upper)
        elif cos_i_fixed is not None:
            self.cos_i = np.array([cos_i_fixed] * self.num_generated)
        else:
            if cos_i_lower is None:
                cos_i_lower = 0.
            if cos_i_upper is None:
                cos_i_upper = 1.
            self.cos_i = np.random.uniform(cos_i_lower, cos_i_upper, self.num_generated)  # uniform distribution between 0 and 1

    def generate_phase(self):
        # Pericenter Phase (radians)
        phase_fixed = self.limits[18]
        phase_lower = self.limits[19]
        phase_upper = self.limits[20]

        # if phase_lower is None and phase_upper is None and phase_fixed is None:
        #     self.phase = np.random.uniform(0, 2*np.pi, self.num_generated)
        if phase_fixed is not None:
            self.phase = np.array([phase_fixed] * self.num_generated)
        else:
            if phase_lower is None:
                phase_lower = 0.
            if phase_upper is None:
                phase_upper = 2. * np.pi
            self.phase = np.random.uniform(phase_lower, phase_upper, self.num_generated)

    def generate_eccentricity(self):
        # Eccentricity
        e_fixed = self.limits[6]
        e_lower = self.limits[7]
        e_upper = self.limits[8]

        log_P = np.log10(self.P)

        if e_fixed is not None:
            self.ecc = np.array([e_fixed] * self.num_generated)
        else:
            # Choose e from within given limits
            a, b, c, d = [0.148, 0.001, 0.042, 0.128]  # parameters from fitting
            dx = 1e-4
            if e_lower is None or e_lower < dx:
                e_lower = dx
            if e_upper is None or e_upper > 0.9999:
                e_upper = 0.9999

            ecc = np.zeros(self.num_generated)
            for i in range(self.num_generated):
                x = np.arange(e_lower, e_upper + dx, dx)  # range of eccentricities
                if log_P[i] < 3:
                    # Calculate the mean and std deviation at this period
                    e_mu = a * log_P[i] + b
                    e_sig = c * log_P[i] + d
                    # Construct the probability distribution
                    #    probability distribution is Gaussian between 0 and 1, and zero outside those bounds
                    y = self.gauss(x, e_mu, e_sig)
                    y = np.divide(y, np.sum(y))
                    # Generate eccentricity
                    ecc[i] = np.random.choice(x, p=y)
                else:
                    ecc[i] = np.random.uniform(e_lower, e_upper)
            self.ecc = ecc

    def generate_arg_peri(self):

        # Argument of Periapsis (radians)
        arg_fixed = self.limits[9]
        arg_lower = self.limits[10]
        arg_upper = self.limits[11]

        if arg_fixed is not None:
            self.arg_peri = np.array([arg_fixed] * self.num_generated)
        else:
            if arg_lower is None:
                arg_lower = 0.
            if arg_upper is None:
                arg_upper = np.pi
            self.arg_peri = np.random.uniform(arg_lower, arg_upper, self.num_generated)

    # Distribution functions
    def gauss(self, x, *p):
        mu, sigma = p
        return np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

    # Accessor function
    def get_all(self):
        return np.vstack([self.P,self.mass_ratio,self.cos_i,self.a,self.ecc,self.arg_peri,self.phase])


    def write(self, filename):
        """ Write companions to an output file
        """
        if (filename.endswith(".hdf5") or filename.endswith(".h5")) == False:
            raise ValueError("This function only supports HDF5 files (.hdf5 or .h5)")

        # By default, want to write to an h5py file. 
        # Applications.py already writes to a .csv (move here eventually)

        nlims = len(self.limits)
        nparams = nlims // 3

        # Open a h5py file ("with")
        with h5py.File(filename,"w") as f:
            # Create a dataset "limits" containing all the prior parameters
            dlims = f.create_dataset("limits",(nlims,),dtype="float32")

            # Populate the limits dataset
            dlims[:] = self.limits[:]

            # Create a dataset with 

            # Create a prior dataset with the existing all_prior values
            all_prior = self.get_all()
            dprior = f.create_dataset("companions",
                                      data=all_prior)



    @classmethod
    def read(cls, filename):
        """ Read in companions from an existing file
        """
        pass