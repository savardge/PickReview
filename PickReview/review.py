import glob
import os
import pyasdf
import obspy
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
from PickReview.utils import *
from PickReview.nlloc_utils import run_nonlinloc
from PickReview.plotting import plot_scatter
import logging
Logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")


def run_snuffler(stream, events, markers, inventory, ntracks=12):
    """
    Launch Pyrocko Snuffler to review picks
    :param stream: Obspy.Stream
    :param events: Obspy Event or Catalog
    :param markers: list of Pyrocko Markers
    :param inventory: Obspy Inventory
    :return: return_tag, list of reviewed Pyrocko Markers
    """
    if isinstance(events, obspy.core.event.Event):
        catalog = obspy.Catalog(events=[event])
    else:
        catalog = events
    if inventory:
        return_tag, markers_out = stream.snuffle(catalog=catalog, inventory=inventory, markers=markers, ntracks=ntracks)
    else:
        return_tag, markers_out = stream.snuffle(catalog=catalog, markers=markers, ntracks=ntracks)
    return return_tag, markers_out


if __name__ == "__main__":

    flist = glob.glob("/home/genevieve/research/quakemigrate/cami/dets_fromscamp_h5/2020031*.h5")
    for filename in flist:

        #print(f"H5 file: {filename}")

        # Open H5 file and extract event and stream
        ds = pyasdf.ASDFDataSet(filename=filename)
        event = ds.events[0].copy()
        stream = obspy.Stream()
        for id in ds.waveforms.list():
            stream += ds.waveforms[id]["raw_recording"]
        del ds

        # Criterion for analysis
        reflat = 50.4502915
        reflon = -112.120833
        refdepth = -543
        tup = gps2dist_azimuth(event.preferred_origin().latitude, event.preferred_origin().longitude, reflat, reflon)
        distance = np.sqrt(tup[0] ** 2 + (event.preferred_origin().depth - refdepth) ** 2)
        if distance > 50:
            continue
        print(event.preferred_origin())

        # Get Pyrocko markers from previous picks
        picks, _ = fix_picks_ids(event, stream, method="modelled")
        markers = picks2markers(picks, event=event, phase=True, kinds=(1, 2))

        # Re-pick using snuffler
        return_tag, markers_out = run_snuffler(stream, events=event, markers=markers, inventory=None)
        #print(f"Return tag: {return_tag}")
        if not return_tag:
            continue

        # Run NonLinLoc with new picks
        run_dir = "/home/genevieve/research/PickReview/NonLinLoc/temp_dir"
        for f in glob.glob(os.path.join(run_dir, "*")):  # Clean out dir
            os.remove(f)
        fname = os.path.join(run_dir, "phases.obs")
        write_obs(phase_markers=markers_out, fname=fname)
        run_nonlinloc(fname)

        # Plot new location
        previous_origin = event.preferred_origin()
        previous_origin.depth += 774
        plot_scatter(hypfile="/home/genevieve/research/PickReview/NonLinLoc/temp_dir/last.hyp")

        # Save new location and picks in a directory (.hyp file)
        savedir = "/home/genevieve/research/PickReview/NonLinLoc/save_dir/"
        dum = os.path.split(glob.glob(os.path.join(run_dir, "result.2*.grid0.loc.hyp"))[0])[1].split(".")
        code = ".".join(dum[1:3])
        os.rename(fname, os.path.join(savedir, "result." + code + ".obs"))

        outfiles = glob.glob(os.path.join(run_dir, "result.*.grid0.loc*"))
        for f in outfiles:
            fc = os.path.split(f)[1]
            if "sum" in fc:
                os.rename(f, os.path.join(savedir, fc.replace("sum", code)))
            else:
                os.rename(f, os.path.join(savedir, fc))