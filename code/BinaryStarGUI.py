# MOLUSC v.20220321
# Mackenna Wood, UNC Chapel Hill
from datetime import datetime
import scipy as scipy
import sys
import warnings
from application import Application
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

today = datetime.today().isoformat().split("T")[0]

if __name__ == '__main__':
    app = Application(sys.argv)
    app.start()
