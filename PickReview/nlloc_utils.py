import glob
import os
from subprocess import Popen, PIPE
import logging
Logger = logging.getLogger(__name__)

NLL_EXECUTABLE = "/home/genevieve/research/PickReview/NonLinLoc/bin/NLLoc"


def run_nonlinloc(fname):
    old_wd = os.getcwd()
    os.chdir(os.path.split(fname)[0])  # Change to NLL run directory
    env = dict(os.environ)

    try:
        p = Popen([NLL_EXECUTABLE + " nlloc_loc.in"], env=env, stdout=PIPE)
    except OSError as e:
        import errno
        if e.errno == errno.ENOENT:
            logging.error('NLLoc executable not found')
            # return
        else:
            raise e

    (out, err) = p.communicate()
    os.chdir(old_wd)

    # Get the output file


#def write_nlloc_control(fout, ):
