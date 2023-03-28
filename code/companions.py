from datetime import datetime
import numpy as np
import scipy as scipy
import scipy.stats as stats
import warnings
from astropy.utils.exceptions import AstropyWarning
import logging
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

        np.random.seed()
        # Period (days), Mass Ratio and Semi-Major Axis (AU)
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

        if a_fixed is None and a_lower is None and a_upper is None:
            # Generate Period independently of mass ratio and separation
            if P_fixed is not None:
                self.P = np.array([P_fixed] * self.num_generated) #TODO: Delete self.P
            else:
                if P_upper is None:
                    P_upper = float('inf')
                self.P = np.array([-1.0] * self.num_generated)  # Initializing
                for i in range(0, self.num_generated):
                    while self.P[i] < P_lower or self.P[i] > P_upper:
                        log_P = np.random.normal(self.mu_log_P, self.sig_log_P) #TODO: Check log_P size for possible deletion
                        self.P[i] = 10 ** log_P            
            
                print("\n\n\nHere are some tables:")
                print(f"log_P: {np.shape(log_P)}\n")
            
            # Now Generate mass ratio
            if mass_fixed is not None:
                self.mass_ratio = np.array([mass_fixed] * self.num_generated) #TODO: Delete self.mass_ratio size for possible deletion
            elif self.q_exp == 0.0:
                # This generates all the  way to q=0 for a uniform distribution
                if mass_lower is None:
                    mass_lower = 0.
                if mass_upper is None:
                    mass_upper = 1.
                self.mass_ratio = np.random.uniform(mass_lower, mass_upper, self.num_generated)
            else:
                # This applies a lower mass limit for the power law distribution
                if mass_lower is None:
                    mass_lower = 0.01
                if mass_upper is None:
                    mass_upper = 1.
                q_range = np.arange(mass_lower, mass_upper + dq, dq)  # range of possible q
                p = np.power(q_range, self.q_exp)    # probabilities
                p = [x / np.sum(p) for x in p]  # normalized probabilities
                self.mass_ratio = np.random.choice(q_range, p=p, size=self.num_generated)
                print("\n\n\nHere are some tables:")
                print(f"q_range: {np.shape(q_range)}\n")
                print(f"p: {np.shape(p)}\n")

            # Now calculate Semi-Major Axis using P and m
            self.a = ((self.P / 365)**2 * G * self.star_mass*(1 + self.mass_ratio)/(4 * np.pi ** 2))**(1/3)  # AU
            # Simple Cases Done
            
            print("\n\n\nHere are some tables:")
            print(f"self.P: {np.shape(self.P)}\n")


            
        elif a_fixed is not None and mass_fixed is not None and P_fixed is not None:
            # Ideally the GUI will not allow you to try to fix all three of them, but for now I will put it here
            return -11
        elif a_fixed is not None:
            if mass_fixed is not None:
                # Set values of a and mass to the given values, calculate P
                self.a = np.array([a_fixed] * self.num_generated)
                self.mass_ratio = np.array([mass_fixed] * self.num_generated)
                self.P = np.sqrt( (4 * np.pi ** 2 * self.a ** 3) / (G * self.star_mass * (1 + self.mass_ratio))) * 365
            elif P_fixed is not None:
                # Set values of a and P to the given values, calculate mass
                self.a = np.array([a_fixed] * self.num_generated)
                self.P = np.array([P_fixed] * self.num_generated)
                self.mass_ratio = np.array(4*np.pi**2 * self.a**3 / (G * self.star_mass * (self.P / 365)** 2) - 1)
                if self.mass_ratio[0] > 1:
                    # The fixed values of period and semi-major axis require a mass ratio greater than one
                    return -12
            else:
                # Only a is fixed.
                self.a = np.array([a_fixed] * self.num_generated)
                if P_lower is None and P_upper is None:
                    # Generate mass ratio as normal and then calculate P
                    # Choose to generate mass ratio instead of P, to ensure that it stays within limits
                    if mass_lower is None:
                        mass_lower = 0.01
                    if mass_upper is None:
                        mass_upper = 1.
                    q_range = np.arange(mass_lower, mass_upper + dq, dq)  # range of possible q
                    p = np.power(q_range, self.q_exp)  # probabilities
                    p = [x / np.sum(p) for x in p]  # normalized probabilities
                    self.mass_ratio = np.random.choice(q_range, p=p, size=self.num_generated)
                    self.P = np.sqrt(
                        (4 * np.pi ** 2 * self.a ** 3) / (G * self.star_mass * (1 + self.mass_ratio))) * 365
                else:
                    # Generate P within limits, and then calculate mass ratio
                    if P_lower is None:
                        P_lower = 0
                    if P_upper is None:
                        P_upper = float('inf')
                    self.P = np.array([-1.0] * self.num_generated)  # Initializing
                    for i in range(0, self.num_generated):
                        while self.P[i] < P_lower or self.P[i] > P_upper:
                            log_P = np.random.normal(self.mu_log_P, self.sig_log_P)
                            self.P[i] = 10 ** log_P
                    self.mass_ratio = np.array(4 * np.pi ** 2 * self.a ** 3 / (G * self.star_mass * (self.P / 365) ** 2) - 1)
        elif a_lower is not None or a_upper is not None:
            if P_fixed is not None:
                self.P = np.array([P_fixed] * self.num_generated)
                # Generate limits on mass based on limits on a
                if mass_lower is None:
                    mass_lower = 0.
                if mass_upper is None:
                    mass_upper = 1.
                # The limits will now need to be a list because the mass ratio for a given separation depends on the P
                mass_lower_list = np.array([0.] * self.num_generated)  # Initialization
                mass_upper_list = np.array([float('inf')] * self.num_generated)  # Initialization

                if a_lower is not None:
                    mass_lower_list = np.array( 4 * np.pi ** 2 * a_lower ** 3 / (G * self.star_mass * (self.P / 365) ** 2) - 1)
                if a_upper is not None:
                    mass_upper_list = np.array(
                        4 * np.pi ** 2 * a_upper ** 3 / (G * self.star_mass * (self.P / 365) ** 2) - 1)

                mass_lower_list = [max(x, mass_lower) for x in mass_lower_list]
                mass_upper_list = [min(x, mass_upper) for x in mass_upper_list]

                self.mass_ratio = np.array([0.] * self.num_generated)
                for i in range(0, self.num_generated):
                    q_range = np.arange(mass_lower, 1 + dq, dq)  # range of possible q
                    p = np.power(q_range, self.q_exp)  # probabilities
                    p = [x / np.sum(p) for x in p]  # normalized probabilities
                    self.mass_ratio[i] = np.random.choice(q_range, p=p, size=1)
                # Calculate a
                self.a = ((self.P / 365) ** 2 * G * self.star_mass * (1 + self.mass_ratio) / (
                            4 * np.pi ** 2)) ** (1 / 3)  # AU
            else:
                # Generate Mass
                if mass_fixed is not None:
                    self.mass_ratio = np.array([mass_fixed] * self.num_generated)
                elif self.q_exp == 0.0:
                    # This generates all the  way to q=0 for a uniform distribution
                    if mass_lower is None:
                        mass_lower = 0.
                    if mass_upper is None:
                        mass_upper = 1.
                    self.mass_ratio = np.random.uniform(mass_lower, mass_upper, self.num_generated)
                else:
                    # This applies a lower mass limit for the power law distribution
                    if mass_lower is None:
                        mass_lower = 0.01
                    if mass_upper is None:
                        mass_upper = 1.
                    q_range = np.arange(mass_lower, mass_upper + dq, dq)  # range of possible q
                    p = np.power(q_range, self.q_exp)  # probabilities
                    p = [x / np.sum(p) for x in p]  # normalized probabilities
                    self.mass_ratio = np.random.choice(q_range, p=p, size=self.num_generated)
                # Generate limits on period, based on limits on a
                if P_lower is None:
                    P_lower = 0.
                if P_upper is None:
                    P_upper = float('inf')

                # The limits will now need to be a list because the period for a given separation depends on the mass
                P_lower_list = np.array([0.] * self.num_generated)  # Initialization
                P_upper_list = np.array([float('inf')] * self.num_generated)  # Initialization

                if a_lower is not None:
                    P_lower_list = np.sqrt(
                        (4 * np.pi ** 2 * a_lower ** 3) / (G * self.star_mass * (1 + self.mass_ratio))) * 365
                if a_upper is not None:
                    P_upper_list = np.sqrt(
                        (4 * np.pi ** 2 * a_upper ** 3) / (G * self.star_mass * (1 + self.mass_ratio))) * 365

                P_lower_list = [max(x, P_lower) for x in P_lower_list]
                P_upper_list = [min(x, P_upper) for x in P_upper_list]

                # Generate period between limits
                self.P = np.array([-1.0] * self.num_generated)  # Initializing
                for i in range(0, self.num_generated):
                    if P_lower_list[i] > P_upper_list[i]:
                        return -11
                    while self.P[i] < P_lower_list[i] or self.P[i] > P_upper_list[i]:
                        log_P = np.random.normal(self.mu_log_P, self.sig_log_P)
                        self.P[i] = 10 ** log_P
                # Calculate a
                self.a = ((self.P / 365)**2 * G * self.star_mass * (1 + self.mass_ratio) / (
                            4 * np.pi ** 2))**(1/3)  # AU
                
                # del(self.mass_ratio) #TODO: Check is deleting self.mass_ratio causes issues with the get_all accessor
        else:
            print("I really don't think you should be able to get here ever...")

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
            #TODO: Check self.cos_i for deletion
            
        print("\n\n\nHere are some tables:")
        print(f"self.cos_i: {np.shape(self.cos_i)}\n")
          
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

        # Eccentricity
        e_fixed = self.limits[6]
        e_lower = self.limits[7]
        e_upper = self.limits[8]

        log_P = np.log10(self.P)

        print("\n\n\nHere are some tables:")
        print(f"log_P: {np.shape(log_P)}\n")

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

        print("\n\n\nHere are some tables:")
        print(f"self.ecc: {np.shape(self.ecc)}\n")
        print(f"x: {np.shape(x)}\n")

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

        print("\n\n\nHere are some tables:")
        print(f"self.arg_peri: {np.shape(self.arg_peri)}\n")
        print(f"self.limits: {np.shape(self.limits)}\n")
        print(f"self.star_mass: {np.shape(self.star_mass)}\n")
        print(f"self.mu_log_P: {np.shape(self.mu_log_P)}\n")
        print(f"self.sig_log_P: {np.shape(self.sig_log_P)}\n")
        print(f"self.q_exp: {np.shape(self.q_exp)}\n\n\n")
        #TODO: Check how large some of these are to see if they are worth deleting after use.

        # del(self.P)


        return

    # Distribution functions
    def gauss(self, x, *p):
        mu, sigma = p
        return np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

    # Accessor function
    def get_all(self):
        return np.vstack([self.P,self.mass_ratio,self.cos_i,self.a,self.ecc,self.arg_peri,self.phase])