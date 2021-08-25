#from utils import *
#from plotting import *
# def run_snuffler(stream, catalog, picks, inv):
#     # Convert picks to markers
#     markers = picks2markers(picks)
#     return_tag, markers_out = stream.snuffle(catalog=catalog, inventory=inv, markers=markers, ntracks=12)
#     return return_tag, markers_out


import glob
import pyasdf
import obspy
from PickReview.utils import *
import logging
Logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")

# flist = glob.glob("/Volumes/ExtremeSSD/quakemigrate/cami/dets_fromscamp_h5/20200310*.h5")
flist = glob.glob("/media/genevieve/ExtremeSSD/quakemigrate/cami/dets_fromscamp_h5/20200310*.h5")
filename = flist[0]
ds = pyasdf.ASDFDataSet(filename=filename)
event = ds.events[0].copy()
stream = obspy.Stream()
for id in ds.waveforms.list():
    stream += ds.waveforms[id]["raw_recording"]
del ds

picks, _ = fix_picks_ids(event, stream, method="modelled")
markers = picks2markers(picks, event=event, phase=True, kinds=(1, 2))

catalog = obspy.Catalog(events=[event])
# return_tag, markers_out = stream.snuffle(catalog=catalog, ntracks=12)
return_tag, markers_out = stream.snuffle(catalog=catalog, markers=markers, ntracks=12)


