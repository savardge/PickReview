import glob
import os
import obspy.core.event
import pyasdf
import obspy
from PickReview.utils import *
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

    # flist = glob.glob("/Volumes/ExtremeSSD/quakemigrate/cami/dets_fromscamp_h5/20200310*.h5")
    flist = glob.glob("/home/genevieve/research/quakemigrate/cami/dets_fromscamp_h5/20200310*.h5")
    filename = flist[0]
    ds = pyasdf.ASDFDataSet(filename=filename)
    event = ds.events[0].copy()
    stream = obspy.Stream()
    for id in ds.waveforms.list():
        stream += ds.waveforms[id]["raw_recording"]
    del ds

    picks, _ = fix_picks_ids(event, stream, method="modelled")
    markers = picks2markers(picks, event=event, phase=True, kinds=(1, 2))

    return_tag, markers_out = run_snuffler(stream, events=event, markers=markers, inventory=None)

    run_dir = "/home/genevieve/research/PickReview/NonLinLoc/temp_dir"
    fname = os.path.join(run_dir, "phases.obs")
    write_obs(phase_markers=markers_out, fname=fname)

