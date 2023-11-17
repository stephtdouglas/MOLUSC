# MOLUSC v.20220321
# Mackenna Wood, UNC Chapel Hill
from datetime import datetime
import os, sys, warnings, logging
import numpy as np

import scipy as scipy
import scipy.stats as stats
from astropy.utils.exceptions import AstropyWarning

# logging.basicConfig(filename=f'{os.getenv("MOLOC")}/../script_logs/molusc.log', format='%(asctime)s %(message)s', encoding='utf-8', level=logging.INFO)

from application import Application

arrayid = int(os.getenv("SLURM_ARRAY_TASK_ID",9999))
jobid = int(os.getenv("SLURM_JOB_ID",9999))
jobname = os.getenv("SLURM_JOB_NAME", 'MOLUSC_999')

logging.basicConfig(level=logging.INFO,
                    filename=f'/data/douglaslab/douglste/script_logs/logger_{jobid}_{arrayid}.log',
                    format='%(asctime)s %(message)s')
logging.getLogger("matplotlib").setLevel(logging.WARNING)


warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.misc.LinAlgWarning)

today = datetime.today().isoformat().split("T")[0]

if __name__ == '__main__':
    app = Application(sys.argv)
    app.start()
