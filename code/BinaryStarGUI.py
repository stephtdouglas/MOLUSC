# MOLUSC v.20220321
# Mackenna Wood, UNC Chapel Hill
from datetime import datetime
import scipy as scipy
import scipy.stats as stats
import sys
import warnings
from application import Application
from astropy.utils.exceptions import AstropyWarning
import logging
import os
# logging.basicConfig(filename=f'{os.getenv("MOLOC")}/../script_logs/molusc.log', format='%(asctime)s %(message)s', encoding='utf-8', level=logging.INFO)

logging.basicConfig(filename=f'molusc.log', format='%(asctime)s %(message)s', encoding='utf-8', level=logging.INFO)
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

today = datetime.today().isoformat().split("T")[0]

if __name__ == '__main__':
    app = Application(sys.argv)
    app.start()