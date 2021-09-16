import glob
import os
import subprocess
import logging
Logger = logging.getLogger(__name__)

NLL_EXECUTABLE = "/home/genevieve/research/PickReview/NonLinLoc/bin/NLLoc"
NLL_CONTROL = "/home/genevieve/research/PickReview/NonLinLoc/run/nlloc_loc.in"


def run_nonlinloc(fname):
    old_wd = os.getcwd()
    os.chdir(os.path.split(fname)[0])  # Change to NLL run directory
    env = dict(os.environ)

    p = subprocess.run([f"{NLL_EXECUTABLE} {NLL_CONTROL}"], env=env, stdout=subprocess.PIPE, shell=True)
    print(p.stdout.decode('utf-8'))
    os.chdir(old_wd)

#def write_nlloc_control(fout, ):
