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

today = datetime.today().isoformat().split("T")[0]

G = 39.478  # Gravitational constant in AU^3/years^2*M_solar

def check_none(x):
    """
    Check if elements of an list or array are None
    Returns an array of booleans
    """

    return np.fromiter((elem is None for elem in x),bool)

class Companions:

    # class functions
    def __init__(self, num_generated, limits, star_mass, pd_mu, pd_sig, q_exp,
                 **params):
        """
        Initialize the Companions class. 

        Required inputs
        ---------------
        num_generated: integer
        limits: list-like, limits/fixed values for parameters. 
            Any of them can be None (note m is q and phi is phase)
            [P_fixed, P_min, P_max, 
            i_fixed, i_min, i_max, e_fixed, e_min, e_max, 
            arg_peri_fixed, arg_peri_min, arg_peri_max,
            m_fixed, m_min, m_max, a_fixed, a_min, a_max, 
            phi_fixed,phi_min, phi_max]
        star_mass: float, primary star mass in solar masses
        pd_mu: median orbital period in log-days
        pd_sig: width of the log-period distribution
        q_exp: power-law slope, if any, for the IMF

        Optional params - these are the orbital parameters read in from a file
        If creating a new Companions object from scratch, do not provide params
        Two of P, q, and a must be provided - the third will then be calculated
        P: array-like, orbital period in days
        q: array-like, mass-ratio
        a: array-like, semi-major axis in AU
        cos_i: array-like, cosine of inclination
        ecc: array-like, eccentricity
        arg_peri: array-like, argument of periastron (small omega)
        phase: array-like, orbital phase (phi)

        """
        self.num_generated = num_generated
        self.limits = limits
        self.star_mass = star_mass
        self.mu_log_P = pd_mu
        self.sig_log_P = pd_sig
        self.q_exp = q_exp

        if params:
            pkeys = params.keys()
            # TODO: this should also take into account fixed values in limits
            if ("P" in pkeys) and ("q" in pkeys):
                self.P = params["P"][:num_generated]
                self.mass_ratio = params["q"][:num_generated]
                if "a" in pkeys:
                    self.a = params["a"][:num_generated]
                else:
                    self.a = ((self.P / 365)**2 * G * self.star_mass*(1 + self.mass_ratio)/(4 * np.pi ** 2))**(1/3)  # AU
            elif ("P" in pkeys) and ("a" in pkeys):
                self.P = params["P"][:num_generated]
                self.a = params["a"][:num_generated]
                if "q" in pkeys:
                    self.mass_ratio = params["q"][:num_generated]
                else:
                    self.mass_ratio = np.array(4*np.pi**2 * self.a**3 / (G * self.star_mass * (self.P / 365)** 2) - 1)
            elif ("a" in pkeys) and ("q" in pkeys):
                self.a = params["a"][:num_generated]
                self.mass_ratio = params["q"][:num_generated]
                if "P" in pkeys:
                    self.P = params["P"][:num_generated]
                else:
                    self.P = np.sqrt( (4 * np.pi ** 2 * self.a ** 3) / (G * self.star_mass * (1 + self.mass_ratio))) * 365
            else:
                raise ValueError("At least two of P, q, and a must be provided")

            # Read in the other parameters
            self.cos_i = params["cos_i"][:num_generated]
            self.ecc = params["ecc"][:num_generated]
            self.arg_peri = params["arg_peri"][:num_generated]
            self.phase = params["phase"][:num_generated]
        else:
            self.P, self.a, self.q = None, None, None
            self.cos_i, self.ecc = None, None 
            self.arg_peri, self.phase = None, None

    def generate(self):
        """
        generate orbital parameters

        """

        np.random.seed(999)

        # If the user has fixed P, a, AND mass_ratio, then we can't continue
        if ((self.limits[0] is not None) and (self.limits[15] is not None) and
            (self.limits[12] is not None)):
            raise ValueError("You can't fix P, a, and q at the same time")

        # If user has not fixed a, then we can generate P and q independently
        elif self.limits[15] is None:
            self.generate_pq()
            self.a = ((self.P / 365)**2 * G * self.star_mass*(1 + self.mass_ratio)/(4 * np.pi ** 2))**(1/3)  # AU

        # If user has fixed a, then P or q becomes a dependent variable
        else:
            self.generate_paq()

        # Then generate all the other orbital parameters
        self.generate_inclination()
        self.generate_phase()
        self.generate_eccentricity()
        self.generate_arg_peri()

    def generate_paq(self):
        """
        generate Period (days), Mass Ratio and Semi-Major Axis (AU)

        """
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


    def generate_pq(self):
        """
        Generate Period (days) and Mass Ratio but NOT semimajor axis

        This function should only be run if 

        """
        a_fixed = self.limits[15]
        if a_fixed is not None:
            raise ValueError("Attempting to run generate_pq with fixed a")

        P_fixed = self.limits[0]
        P_lower = self.limits[1]
        P_upper = self.limits[2]
        mass_fixed = self.limits[12]
        mass_lower = self.limits[13]
        mass_upper = self.limits[14]

        print('Generation Mass Ratio Exponent: ', self.q_exp)
        dq = 1e-4  # coarseness

        # Default lower limit on P
        if P_lower is None or P_lower < 0.1:
            P_lower = 0.1

        # We will always generate P. Either it's fixed, OR the
        # other two values are fixed, OR we generate it independently
        if P_fixed is not None:
            self.P = np.full(self.num_generated,fill_value=P_fixed)

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

        # We've handled all cases where values can be fixed. 
        # Now we generate the mass ratio randomly if needed, 
        # and calculate a if it wasn't fixed
        if mass_fixed is not None:
            self.mass_ratio = np.full(self.num_generated,fill_value=mass_fixed)

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
            # Create a group to contain metadata and misc parameters
            drp = f.create_group("meta")

            # Create a dataset "limits" containing all the prior parameters
            dlims = drp.create_dataset("limits",(nlims,),dtype="float32")

            # Populate the limits dataset
            dlims[:] = self.limits[:]

            # Add the other misc parameters
            drp.create_dataset("star_mass",data=self.star_mass)
            drp.create_dataset("mu_log_P",data=self.mu_log_P)
            drp.create_dataset("sig_log_P",data=self.sig_log_P)
            drp.create_dataset("q_exp",data=self.q_exp)
            drp.create_dataset("num_generated",data=self.num_generated)

            # Create a group to contain the existing prior values

            grp = f.create_group("companions")

            # Write out the parameters by name
            P_lims = self.limits[:2]
            mass_lims = self.limits[12:15]
            a_lims = self.limits[15:18]
            P_check = np.any(check_none(P_lims)==False)
            mass_check = np.any(check_none(mass_lims)==False)
            a_check = np.any(check_none(a_lims)==False)
            # If 2 values were fixed, that determines the other by default
            # Just ignore the unfixed one and it will be re-calculated
            # The only case where P was separately calculated was when a and q were fixed
            if a_check and mass_check:
                grp.create_dataset("q",data=self.mass_ratio)
                grp.create_dataset("a",data=self.a)

            # If only a was fixed, then P was drawn randomly and the 
            # stellar mass informs the mass ratio
            elif a_check:
                grp.create_dataset("P",data=self.P)
                grp.create_dataset("a",data=self.a)


            # In any other cases, we calculate P and q randomly or one
            # of them is fixed, but either way they're the only two we need
            else:
                grp.create_dataset("P",data=self.P)
                grp.create_dataset("q",data=self.mass_ratio)


            grp.create_dataset("cos_i",data=self.cos_i)
            grp.create_dataset("ecc",data=self.ecc)
            grp.create_dataset("arg_peri",data=self.arg_peri)
            grp.create_dataset("phase",data=self.phase)



    @classmethod
    def read(cls, filename, star_mass=None, num_max=None, warn_over=True):
        """ Read in companions from an existing file and initialize
        the companions object.

        Things that can change: target mass, max number of companions

        Things that definitely can't change: mu/sigma for the logP distribution
            slope of the IMF, mass_ratio, cos_i, ecc, arg_peri, phase

        ???????: P, a
        """

        if (filename.endswith(".hdf5") or filename.endswith(".h5")) == False:
            raise ValueError("This function only supports HDF5 files (.hdf5 or .h5)")

        with h5py.File(filename,"r") as f:

            ngen = np.int64(f["meta"]["num_generated"][()])

            if num_max is None:
                nread = ngen
            else:
                if num_max <= ngen:
                    nread = num_max
                elif warn_over and (num_max>ngen):
                    raise ValueError("Requested more companions than are provided in the file")
                else:
                    nread = ngen

            # If a new mass is not provided, assume we're using the one
            # from the hdf5 file
            if star_mass is None:
                star_mass = f["meta"]["star_mass"][()]


            return cls(nread,f["meta"]["limits"][:],
                       star_mass,f["meta"]["mu_log_P"][()],
                       f["meta"]["sig_log_P"][()],f["meta"]["q_exp"][()],
                       **f["companions"])

