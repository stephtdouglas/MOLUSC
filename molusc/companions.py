from datetime import datetime
from collections import OrderedDict
import warnings
import logging

import h5py
import numpy as np
import scipy as scipy
import scipy.stats as stats
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.LinAlgWarning)

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
        This is the prior for rejection sampling.

        Required inputs
        ---------------
        num_generated: integer
        limits: dictionary, limits/fixed values for parameters. 
            dictionary keys are parameters, 
            values are OrderedDict of "fixed","min","max"
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
                    self._calc_a_from_Pq()
            elif ("P" in pkeys) and ("a" in pkeys):
                self.P = params["P"][:num_generated]
                self.a = params["a"][:num_generated]
                if "q" in pkeys:
                    self.mass_ratio = params["q"][:num_generated]
                else:
                    self._calc_q_from_aP()
            elif ("a" in pkeys) and ("q" in pkeys):
                self.a = params["a"][:num_generated]
                self.mass_ratio = params["q"][:num_generated]
                if "P" in pkeys:
                    self.P = params["P"][:num_generated]
                else:
                    self._calc_P_from_aq()
            else:
                raise ValueError("At least two of P, q, and a must be provided")

            # Read in the other parameters
            self.cos_i = params["cos_i"][:num_generated]
            self.ecc = params["ecc"][:num_generated]
            self.arg_peri = params["arg_peri"][:num_generated]
            self.phase = params["phase"][:num_generated]

            if "v0" in pkeys:
                self.v0 = params["v0"][:num_generated]
            else:
                self.v0 = None
        else:
            self.P, self.a, self.q = None, None, None
            self.cos_i, self.ecc = None, None 
            self.arg_peri, self.phase = None, None
            self.v0 = None

    def generate(self):
        """
        generate orbital parameters

        """

        np.random.seed(999)

        # If the user has fixed P, a, AND mass_ratio, then we can't continue
        P_lims = np.fromiter(self.limits["P"].values(),"float")
        a_lims = np.fromiter(self.limits["a"].values(),"float")
        q_lims = np.fromiter(self.limits["q"].values(),"float")
        # These check params will be True if there are any non-null values
        # If no limits/fixed values are set, then the check params are False
        P_check = np.any(np.isfinite(P_lims))
        a_check = np.any(np.isfinite(a_lims))
        q_check = np.any(np.isfinite(q_lims))
        if (P_check and a_check and q_check):
            raise ValueError("You can't fix P, a, and q at the same time")

        # If user has not fixed a, then we can generate P and q independently
        elif a_check==False:
            self.generate_pq()
            self._calc_a_from_Pq()

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

        print('Generation Mass Ratio Exponent: ', self.q_exp)
        dq = 1e-4  # coarseness

        # Default lower limit on P
        if self.limits["P"]["min"] is None or self.limits["P"]["min"] < 0.1:
            P_lower = 0.1
        else:
            P_lower = self.limits["P"]["min"]

        if self.limits["q"]["fixed"] is not None:
            self.mass_ratio = np.full(self.num_generated,fill_value=self.limits["q"]["fixed"])

        if self.limits["a"]["fixed"] is not None:
            self.a = np.full(self.num_generated,fill_value=self.limits["a"]["fixed"])

        # We will always generate P. Either it's fixed, OR the
        # other two values are fixed, OR we generate it independently
        if self.limits["P"]["fixed"] is not None:
            self.P = np.full(self.num_generated,
                             fill_value=self.limits["P"]["fixed"])
        # with a and q, calculate P directly
        elif ((self.limits["a"]["fixed"] is not None) and 
              (self.limits["q"]["fixed"] is not None)):
            self._calc_P_from_aq()
            
        else:
            tn_low = (P_lower-self.mu_log_P)/self.sig_log_P
            if self.limits["P"]["max"] is None:
                tn_up = np.inf
            else:
                tn_up = (self.limits["P"]["max"]-self.mu_log_P)/self.sig_log_P
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
        if (self.limits["a"]["fixed"] is not None):
            self._calc_q_from_aP()
            if np.any(self.mass_ratio > 1):
                # The fixed values of period and semi-major axis require a mass ratio greater than one
                return -12

        # We've handled all cases where values can be fixed. 
        # Now we generate the mass ratio randomly if needed, 
        # and calculate a if it wasn't fixed
        else:
            # We already handled this case above
            if self.limits["q"]["fixed"] is not None:
                pass

            # Flat IMF
            elif self.q_exp == 0.0:
                # This generates all the  way to q=0 for a uniform distribution
                if self.limits["q"]["min"] is None:
                    mass_lower = 0.
                else:
                    mass_lower = self.limits["q"]["min"]
                if self.limits["q"]["max"] is None:
                    mass_upper = 1.
                else:
                    mass_upper = self.limits["q"]["max"]
                self.mass_ratio = np.random.uniform(mass_lower, mass_upper, self.num_generated)

            # IMF with an input slope
            else:
                # This applies a lower mass limit for the power law distribution
                if self.limits["q"]["min"] is None:
                    mass_lower = 0.01
                else:
                    mass_lower = self.limits["q"]["min"]
                if self.limits["q"]["max"] is None:
                    mass_upper = 1.
                else:
                    mass_upper = self.limits["q"]["max"]
                q_range = np.arange(mass_lower, mass_upper + dq, dq)  # range of possible q
                p = np.power(q_range, self.q_exp)    # probabilities
                p = p/np.sum(p) # normalized probabilities
                self.mass_ratio = np.random.choice(q_range, p=p, size=self.num_generated)
            # Now calculate Semi-Major Axis using P and q
            self._calc_a_from_Pq()       


    def generate_pq(self):
        """
        Generate Period (days) and Mass Ratio but NOT semimajor axis

        This function should only be run if a is not fixed in any way

        """
        a_fixed = self.limits["a"]["fixed"]
        if a_fixed is not None:
            raise ValueError("Attempting to run generate_pq with fixed a")

        print('Generation Mass Ratio Exponent: ', self.q_exp)
        dq = 1e-4  # coarseness

        # Default lower limit on P
        if self.limits["P"]["min"] is None or self.limits["P"]["min"] < 0.1:
            P_lower = 0.1
        else:
            P_lower = self.limits["P"]["min"]

        # We will always generate P. Either it's fixed, OR the
        # other two values are fixed, OR we generate it independently
        if self.limits["P"]["fixed"] is not None:
            self.P = np.full(self.num_generated,
                             fill_value=self.limits["P"]["fixed"])

        else:
            tn_low = (P_lower-self.mu_log_P)/self.sig_log_P
            if self.limits["P"]["max"] is None:
                tn_up = np.inf
            else:
                tn_up = (elf.limits["P"]["max"]-self.mu_log_P)/self.sig_log_P
            P_tn = stats.truncnorm(tn_low,tn_up,
                                   loc=self.mu_log_P,scale=self.sig_log_P)
            log_P = P_tn.rvs(self.num_generated)
            self.P = 10**log_P

        # We've handled all cases where values can be fixed. 
        # Now we generate the mass ratio randomly if needed, 
        # and calculate a if it wasn't fixed
        if self.limits["q"]["fixed"] is not None:
            self.mass_ratio = np.full(self.num_generated,
                                      fill_value=self.limits["q"]["fixed"])

        # Flat IMF
        elif self.q_exp == 0.0:
            # This generates all the  way to q=0 for a uniform distribution
            if self.limits["q"]["min"] is None:
                mass_lower = 0.
            else:
                mass_lower = self.limits["q"]["min"]
            if self.limits["q"]["max"] is None:
                mass_upper = 1.
            else:
                mass_upper = self.limits["q"]["max"]
            self.mass_ratio = np.random.uniform(mass_lower, mass_upper, self.num_generated)

        # IMF with an input slope
        else:
            # This applies a lower mass limit for the power law distribution
            if self.limits["q"]["min"] is None:
                mass_lower = 0.01
            else:
                mass_lower = self.limits["q"]["min"]
            if self.limits["q"]["max"] is None:
                mass_upper = 1.
            else:
                mass_upper = self.limits["q"]["max"]
            q_range = np.arange(mass_lower, mass_upper + dq, dq)  # range of possible q
            p = np.power(q_range, self.q_exp)    # probabilities
            p = p/np.sum(p) # normalized probabilities
            self.mass_ratio = np.random.choice(q_range, p=p, size=self.num_generated)
      
            

    def generate_inclination(self):
        # Inclination
        cos_i_fixed = self.limits["cos_i"]["fixed"]
        cos_i_lower = self.limits["cos_i"]["min"]
        cos_i_upper = self.limits["cos_i"]["max"]

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
            self.cos_i = np.full(self.num_generated,fill_value=cos_i_fixed)
        else:
            if cos_i_lower is None:
                cos_i_lower = 0.
            if cos_i_upper is None:
                cos_i_upper = 1.
            self.cos_i = np.random.uniform(cos_i_lower, cos_i_upper, self.num_generated)  # uniform distribution between 0 and 1

    def generate_phase(self):
        # Pericenter Phase (radians)
        phase_fixed = self.limits["phase"]["fixed"]
        phase_lower = self.limits["phase"]["min"]
        phase_upper = self.limits["phase"]["max"]

        # if phase_lower is None and phase_upper is None and phase_fixed is None:
        #     self.phase = np.random.uniform(0, 2*np.pi, self.num_generated)
        if phase_fixed is not None:
            self.phase = np.full(self.num_generated,fill_value=phase_fixed)
        else:
            if phase_lower is None:
                phase_lower = 0.
            if phase_upper is None:
                phase_upper = 2. * np.pi
            self.phase = np.random.uniform(phase_lower, phase_upper, self.num_generated)

    def generate_eccentricity(self):
        # Eccentricity
        e_fixed = self.limits["ecc"]["fixed"]
        e_lower = self.limits["ecc"]["min"]
        e_upper = self.limits["ecc"]["max"]

        log_P = np.log10(self.P)

        if e_fixed is not None:
            self.ecc = np.full(self.num_generated,fill_value=e_fixed)
        else:
            # Choose e from within given limits
            a, b, c, d = [0.148, 0.001, 0.042, 0.128]  # parameters from fitting
            dx = 1e-4
            if e_lower is None or e_lower < dx:
                e_lower = dx
            if e_upper is None or e_upper > 0.9999:
                e_upper = 0.9999

            ecc = np.zeros(self.num_generated)
            x = np.arange(e_lower, e_upper + dx, dx)  # range of eccentricities
            for i in range(self.num_generated):
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
        arg_fixed = self.limits["arg_peri"]["fixed"]
        arg_lower = self.limits["arg_peri"]["min"]
        arg_upper = self.limits["arg_peri"]["max"]

        if arg_fixed is not None:
            self.arg_peri = np.full(self.num_generated,fill_value=arg_fixed)
        else:
            if arg_lower is None:
                arg_lower = 0.
            if arg_upper is None:
                arg_upper = np.pi
            self.arg_peri = np.random.uniform(arg_lower, arg_upper, self.num_generated)

    def generate_v0(self):

        # system velocity in km/s
        v0_fixed = self.limits["v0"]["fixed"]
        # v0_lower = self.limits["v0"]["min"]
        # v0_upper = self.limits["v0"]["max"]
        v0_mu = self.limits["v0"]["mu"]
        v0_sigma = self.limits["v0"]["sigma"]

        if v0_fixed is not None:
            self.v0 = np.full(self.num_generated,fill_value=v0_fixed)
        elif (v0_mu is not None) and (v0_sigma is not None):
            self.v0 = np.random.normal(loc=v0_mu,scale=v0_sigma,
                                       size=self.num_generated)
        else:
            # For now, anything else means we do a fit as before
            self.v0 = None
            

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

        # Open a h5py file ("with")
        with h5py.File(filename,"w") as f:

            # Create a dataset "limits" containing all the prior parameters
            dlims = f.create_group("limits")

            # Populate the limits group
            for key, val in self.limits.items():
                limlist = np.fromiter(val.values(),"float")
                dlims.create_dataset(key,data=limlist)



            # Create a group to contain metadata and misc parameters
            drp = f.create_group("meta")
            drp.create_dataset("star_mass",data=self.star_mass)
            drp.create_dataset("mu_log_P",data=self.mu_log_P)
            drp.create_dataset("sig_log_P",data=self.sig_log_P)
            drp.create_dataset("q_exp",data=self.q_exp)
            drp.create_dataset("num_generated",data=self.num_generated)

            # Create a group to contain the existing prior values

            grp = f.create_group("companions")

            # Write out the parameters by name
            P_lims = np.fromiter(self.limits["P"].values(),"float")
            a_lims = np.fromiter(self.limits["a"].values(),"float")
            q_lims = np.fromiter(self.limits["q"].values(),"float")
            P_check = np.any(np.isfinite(P_lims))
            a_check = np.any(np.isfinite(a_lims))
            q_check = np.any(np.isfinite(q_lims))
            # If 2 values were fixed, that determines the other by default
            # Just ignore the unfixed one and it will be re-calculated
            # The only case where P was separately calculated was when a and q were fixed
            if a_check and q_check:
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


            # Read limits back into a dictionary
            limits = {}
            sub_keys = ["fixed","min","max"]
            dict_group_keys = f["limits"].keys()
            for key in dict_group_keys:
                limits[key] = OrderedDict()
                for i, sub_key in enumerate(sub_keys):
                    if np.isnan(f["limits"][key][i]):
                        limits[key][sub_key] = None
                    else:
                        limits[key][sub_key] = f["limits"][key][i]

            return cls(nread,limits,
                       star_mass,f["meta"]["mu_log_P"][()],
                       f["meta"]["sig_log_P"][()],f["meta"]["q_exp"][()],
                       **f["companions"])

    def _calc_P_from_aq(self):
        """
        Calculate P when a and q are already set
        """
        self.P = np.sqrt( (4 * np.pi ** 2 * self.a ** 3) / (G * self.star_mass * (1 + self.mass_ratio))) * 365

    def _calc_a_from_Pq(self):
        """
        Calculate a when a and P are already set
        """
        self.a = ((self.P / 365)**2 * G * self.star_mass*(1 + self.mass_ratio)/(4 * np.pi ** 2))**(1/3)  # AU

    def _calc_q_from_aP(self):
        """
        Calculate q when a and P are already set
        """
        self.mass_ratio = np.array(4*np.pi**2 * self.a**3 / (G * self.star_mass * (self.P / 365)** 2) - 1)

