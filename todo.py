from __future__ import print_function

import os
import tempfile
import math
import glob
import numpy as num
import shutil
import datetime

from pyrocko.gui.snuffling import Snuffling, Switch, Param, Choice
from pyrocko.gui.pile_viewer import EventMarker, PhaseMarker
from pyrocko import util, model, orthodrome
from pyrocko.plot import gmtpy
from subprocess import Popen, PIPE

from obspy.core.event import Catalog, Event, Pick, CreationInfo, WaveformStreamID
from obspy.core.event.base import QuantityError, Comment
from obspy.core.event.origin import Origin
from obspy import UTCDateTime

deg2rad = math.pi / 180.
pjoin = os.path.join
reftime = datetime.datetime(1970, 1, 1)


def nsl_str(nsl):
    return '.'.join(nsl)


class NonLinLoc(Snuffling):
    '''
    <html>
    <head>
    <style type="text/css">
        body { margin-left:10px };
    </style>
    </head>
    <body>
    <h1 align="center">Locate seismic sources using NonLinLoc</h1>
    Velocity grids and time grids are assumed to be pre-computed for all relevant stations.

    <h2>Usage</h2>
    &middot; Invoke snuffler with station meta data.
    <br>
    &middot; Pick <i>P</i> and <i>S</i> phases and define them as
    such by pressing <i>F1</i> and <i>F2</i>, respectively.
    <br>
    &middot; Select markers you want to use in the inversion (e.g.: press
    <i>a</i> to select all markers in the visible time range)
    <br>
    &middot; Press <i>Run</i>
    <p>
    The ouput is printed to the terminal. For further information on

    </body>
    </html>
    '''

    def setup(self):

        self.set_name('Nonlinloc location')
        self.nonlinloc_data_dir = pjoin(self._path, 'nonlinloc/data')
        self.add_parameter(
            Param('P std. deviation [s]', 'p_stddev', 0.1, 0.001, 0.2))
        self.add_parameter(
            Param('S std. deviation [s]', 's_stddev', 0.2, 0.001, 0.2))
        self.event_obspy = None

        self.setup_gui(reloaded=True)
        self.set_live_update(True)
        self.dir = None

    def call(self):
        '''Main work routine of the snuffling.'''

        self.cleanup()

        viewer = self.get_viewer()

        self.markers = self.get_selected_markers()
        if len(self.markers) == 0:
            self.fail('No markers selected.')

        self.active_event = viewer.get_active_event()

        if len(viewer.stations) == 0:
            self.fail('No station information available.')

        self.markers2obspy()

        # obspy event to NonLinLoc Input
        cat = Catalog()
        cat.append(self.event_obspy)
        self.dir = tempfile.mkdtemp(prefix='nlloc-%s-' % os.environ['USER'])
        fname = pjoin("/Users/genevieve/contrib-snufflings/nonlinloc-snuffling/nonlinloc", 'picks.obs')
        print(fname)
        print(self.event_obspy)
        self.event_obspy.write(fname, format="NLLOC_OBS")

        print()
        print('=== Running NonLinLoc ' + '=' * 80)
        print('temp dir: %s' % self.dir)
        print('obs file name: %s' % fname)

        old_wd = os.getcwd()
        os.chdir("/Users/genevieve/contrib-snufflings/nonlinloc-snuffling/nonlinloc")
        env = dict(os.environ)

        try:
            p = Popen([pjoin(self._path, 'nonlinloc/bin/NLLoc nlloc.in')], env=env, stdout=PIPE)
        except OSError as e:
            try:
                # try included, compiled version:
                abs_path = os.path.dirname(os.path.abspath(__file__))
                executable = pjoin(abs_path, 'nonlinloc', 'bin', 'NLLoc')
                print('abspath', executable)
                p = Popen([executable], env=env, stdout=PIPE)
            except OSError as e:
                import errno
                if e.errno == errno.ENOENT:
                    self.fail('NLLoc executable not found')
                    #return
                else:
                    raise e

        (out, err) = p.communicate()
        os.chdir(old_wd)

        from glob import glob
        flist = glob(pjoin("/Users/genevieve/contrib-snufflings/nonlinloc-snuffling/nonlinloc/picks*.hyp"))
        print(flist)
        if flist and len(flist) == 1:
            hypfile = flist[0]
            self.plot_scatter(hypfile=hypfile)


    def plot_scatter(self, hypfile, picked_stations):
        from obspy.io.nlloc.util import read_nlloc_scatter
        from obspy import read_events
        import numpy as np
        from pyproj import Proj

        figfile = os.path.join(os.path.split(hypfile)[1].split("grid0")[0] + "_map.png")
        cat = read_events(hypfile, format="NLLOC_HYP")
        picked_stations = list(set([p.waveform_id.station_code for p in cat[0].picks]))
        print("Picked stations: ")
        print(picked_stations)

        # lat, long to X,Y
        injwell = (50.45041, -112.12067)
        p = Proj(proj='eqc', lat_0=injwell[0], lon_0=injwell[1])
        evlat = cat[0].origins[0].latitude
        evlon = cat[0].origins[0].longitude
        xev, yev = p(evlon, evlat)
        xev *= 1e-3
        yev *= 1e-3
        zev = cat[0].origins[0].depth * 1e-3
        print(xev, yev, zev)

        scatter_cloud = read_nlloc_scatter(hypfile.replace(".hyp", ".scat"))
        scatter_cloud = sorted(scatter_cloud, key=lambda x: x[3])

        x = [s[0] for s in scatter_cloud]
        y = [s[1] for s in scatter_cloud]
        z = [s[2] for s in scatter_cloud]
        pdf = [s[3] for s in scatter_cloud]

        xbest = x[-1]
        ybest = y[-1]
        zbest = z[-1]

        import pandas as pd
        stalst = pd.read_csv("stations_XYZ.csv")
        xs = stalst["x"].values
        ys = stalst["y"].values
        zs = stalst["z"].values
        xsbh = np.append(stalst.loc[stalst.station.str.startswith("G")].x.values, 0)
        ysbh = np.append(stalst.loc[stalst.station.str.startswith("G")].y.values, 0)
        zsbh = np.append(stalst.loc[stalst.station.str.startswith("G")].z.values, 0)
        xspick = []
        yspick = []
        zspick = []
        for sta in picked_stations:
            xspick.append(stalst.loc[stalst.station == sta].x.values)
            yspick.append(stalst.loc[stalst.station == sta].y.values)
            zspick.append(stalst.loc[stalst.station == sta].z.values)

        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 14})

        #fig = self.pylab(get='figure')
        fig, axs = plt.subplots(2, 2, figsize=(14, 14))
        axs[1][1].set_axis_off()

        # X-Y
        axs[0][0].scatter(x, y, c=pdf)
        axs[0][0].scatter(xbest, ybest, s=150, c="r", marker="*", edgecolors="k", linewidths=1.0)
        # axs[0][0].scatter(xev, yev, s=50, c="c", marker="*")
        axs[0][0].scatter(xs, ys, s=50, c="c", marker="s", edgecolors="k", linewidths=1.0, zorder=40)
        axs[0][0].scatter(xspick, yspick, s=50, c="k", marker=".", zorder=41)
        axs[0][0].plot(0, 0, marker="o", markersize=15, color="m", markeredgecolor="k",
                       markeredgewidth=4)  # Injection well
        axs[0][0].set_xlabel("Easting (km)")
        axs[0][0].set_ylabel("Northing (km)")
        axs[0][0].xaxis.set_label_position('top')

        # X-Z
        # geology
        axs[1][0].axhline(0, -1, 1, color="k", zorder=1)
        axs[1][0].axhspan(0.025, 0.099, facecolor="brown", alpha=0.25, zorder=1, ec="k",
                          lw=1)  # bedrock/Dinosaur Park top
        axs[1][0].axhspan(0.0415, 0.046, facecolor="k", alpha=0.25, zorder=2, ec="k", lw=1)  # coal zone
        axs[1][0].axhspan(0.099, 0.143, facecolor="green", alpha=0.25, zorder=1, ec="k", lw=1)  # Oldman Fm
        axs[1][0].axhspan(0.157, 0.160, facecolor="k", alpha=0.25, zorder=2, ec="k", lw=1)  # coal zone Taber
        axs[1][0].axhspan(0.271, 0.295, facecolor="k", alpha=0.25, zorder=2, ec="k", lw=1)  # coal zone MacKay
        axs[1][0].axhspan(0.295, 0.302, facecolor="m", alpha=0.25, zorder=2, ec="k", lw=1)  # BBRS
        axs[1][0].plot([0, 0], [0, 0.35], "k", linewidth=2)
        # locs
        axs[1][0].scatter(x, z, c=pdf, zorder=20)
        axs[1][0].scatter(xbest, zbest, s=150, c="r", marker="*", edgecolors="k", linewidths=1.0, zorder=30)
        # stations
        axs[1][0].scatter(xs, zs, s=50, c="c", marker="s", edgecolors="k", linewidths=1.0, zorder=40)
        axs[1][0].scatter(xspick, zspick, s=50, c="k", marker=".", zorder=41)
        axs[1][0].plot(0, 0, marker="o", markersize=15, color="m", markeredgecolor="k",
                       markeredgewidth=4)  # Injection well
        axs[1][0].invert_yaxis()
        axs[1][0].set_ylabel("Depth (km)")
        # ax.set_yticklabels([])
        # ax.set_xticklabels([])

        # Y-Z
        # geology
        axs[0][1].axvline(0, -1, 1, color="k", zorder=1)
        axs[0][1].axvspan(0.025, 0.099, facecolor="brown", alpha=0.25, zorder=1, ec="k",
                          lw=1)  # bedrock/Dinosaur Park top
        axs[0][1].axvspan(0.0415, 0.046, facecolor="k", alpha=0.25, zorder=2, ec="k", lw=1)  # coal zone
        axs[0][1].axvspan(0.099, 0.143, facecolor="green", alpha=0.25, zorder=1, ec="k", lw=1)  # Oldman Fm
        axs[0][1].axvspan(0.157, 0.160, facecolor="k", alpha=0.25, zorder=2, ec="k", lw=1)  # coal zone Taber
        axs[0][1].axvspan(0.271, 0.295, facecolor="k", alpha=0.25, zorder=2, ec="k", lw=1)  # coal zone MacKay
        axs[0][1].axvspan(0.295, 0.302, facecolor="m", alpha=0.25, zorder=2, ec="k", lw=1)  # BBRS
        axs[0][1].plot([0, 0.35], [0, 0], "k", linewidth=2)
        # locs
        axs[0][1].scatter(z, y, c=pdf, zorder=20)
        axs[0][1].scatter(zbest, ybest, s=150, c="r", marker="*", edgecolors="k", linewidths=1.0, zorder=30)
        # stations/well
        axs[0][1].scatter(zs, ys, s=50, c="c", marker="s", edgecolors="k", linewidths=1.0, zorder=40)
        axs[0][1].scatter(zspick, yspick, s=50, c="k", marker=".", zorder=41)
        axs[0][1].plot(0, 0, marker="o", markersize=15, color="m", markeredgecolor="k",
                       markeredgewidth=4)  # Injection well
        axs[0][1].set_yticklabels([])
        axs[0][1].set_xlabel("Depth (km)")

        plt.subplots_adjust(wspace=0.01, hspace=0.01)
        plt.show()
        #plt.savefig(figfile)
        #plt.close()

    def markers2obspy(self):
        creation_string = CreationInfo(creation_time=UTCDateTime.now())
        event_time = self.active_event.time
        latlon_error = QuantityError(uncertainty=0.009)  # equals 1 km
        event = Event()

        for marker in self.markers:
            if not isinstance(marker, PhaseMarker) or self.use_active and \
                    marker.get_event() != self.active_event:
                continue
            phasename = marker.get_phasename()
            nslcs = list(marker.nslc_ids)
            station = nslcs[0][1]

            phase_station = list(marker.nslc_ids)[0][1]
            phase_network = list(marker.nslc_ids)[0][0]
            phase_channel = list(marker.nslc_ids)[0][3]
            wav_id = WaveformStreamID(station_code=phase_station,
                                      channel_code=phase_channel,
                                      network_code=phase_network)
            phase_name = marker.get_phasename()
            phase_tmin = marker.tmin
            phase_tmax = marker.tmax
            if phase_tmin == phase_tmax:
                if phasename == 'P':
                    phase_unc = self.p_stddev
                elif phasename == 'S':
                    phase_unc = self.s_stddev
                phase_time = phase_tmin
            else:
                phase_unc = (phase_tmax - phase_tmin) * 0.5
                phase_mid_int = (phase_tmax - phase_tmin) * 0.5
                phase_time = phase_tmin + phase_mid_int

            phase_time_utc = UTCDateTime(reftime.utcfromtimestamp(phase_time))
            if marker.get_polarity() == 1:
                phase_polarity = "positive"
            elif marker.get_polarity() == -1:
                phase_polarity = "negative"
            else:
                phase_polarity = "undecidable"

            pick = Pick(time=phase_time_utc, phase_hint=phase_name, waveform_id=wav_id, polarity=phase_polarity,
                        time_errors=QuantityError(uncertainty=phase_unc), method_id="snuffler")
            event.picks.append(pick)

        self.event_obspy = event

    def save_last_run(self):
        if not self.dir:
            self.fail('Run first')
        outdir = self.output_filename(caption='Save directory')
        shutil.move(self.dir, outdir)


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    return [NonLinLoc()]


# -*- coding: utf-8 -*-

import datetime
import os

from pyrocko.snuffling import Choice, Snuffling, Switch
from pyrocko.pile_viewer import Marker


# Defaults to user's home folder. os.path.expanduser('~')
USER_DEF_ROOT_FOLDER = \
    os.path.expanduser('~')
USER_DEF_PICKING_FOLDER = "/Volumes/SanDisk_2TB/pyrocko_files/snuffler_picks"


class QuickSave(Snuffling):
    """
    <html>
    <h2>Quick Save</h2>
    <body>
    <p>
    Choose the predefined or user defined path from the <b>Root folder</b>
    menu, desired <b>options</b> and press <b>Run</b>. Path to marker file is
    printed to the terminal.<br>
    User defined path is selected by default and can be modified in the source
    file. It is given as a string under <i>USER_DEF_ROOT_FOLDER</i> variable
    (line 11). Its default value points to user's home folder.
    </p>
    <p>
    Options:</br>
    <ul>
    <li>Save only selected markers.</li>
    <li>Save only markers associated to the active event.</li>
        <ul>
        <li>Add active event name to path. Creates new folder if needed.</li>
        <li>Use active event name as filename.</li>
        <li>Save to Velest format.</li>
        </ul>
    </ul>
    </p>
    <p>
    Author: G. Rajh.
    </p>
    </body>
    </html>
    """

    def setup(self):
        self.set_name('Quick Save')

        home_folder = os.path.expanduser('~')

        self.add_parameter(Choice('Root folder', 'root_folder',
                                  USER_DEF_ROOT_FOLDER, ["/Volumes/SanDisk_2TB/eqcorrscan/manual_picks/repick",
                                                         USER_DEF_PICKING_FOLDER,
                                                         home_folder + '/Desktop',
                                                         home_folder + '/cami/pyrocko',
                                                         home_folder + '/frs_borehole/pyrocko',
                                                         home_folder + '/cami/nonlinloc/obs',
                                                         '/Volumes/SanDisk_2TB/picking'
                                                         ]))

        self.add_parameter(Switch(
            'Save only selected markers', 'save_sel_markers', False))
        self.add_parameter(Switch(
            'Save only markers associated\nto the active event',
            'save_asc_markers', False))
        self.add_parameter(Switch('Add active event name to path',
                                  'name_to_path', False))
        self.add_parameter(Switch('Use active event name as filename',
                                  'name_as_filename', False))
        self.add_parameter(Switch('Save to Velest format',
                                  'velest_format', False))
        self.add_parameter(Switch('Save to QuakeML & NLLOC format',
                                  'quakeml_format', False))

        self.setup_gui(reloaded=True)
        self.set_live_update(False)

    def folder_check(self, save_path):
        try:
            os.makedirs(save_path)
        except OSError:
            if not os.path.isdir(save_path):
                raise

    def check_params(self):
        if self.save_asc_markers and self.save_sel_markers:
            self.error(' Select either "Save only selected markers" or \
"Save only markers associated to the active event" option.')
            self.checked = False

        elif not self.save_asc_markers and self.name_to_path:
            self.save_asc_markers = True
            self.save_sel_markers = False
            self.checked = True

        elif not self.save_asc_markers and self.name_as_filename:
            self.save_asc_markers = True
            self.save_sel_markers = False
            self.checked = True

        elif not self.save_asc_markers and (self.velest_format or self.quakeml_format):
            self.save_asc_markers = True
            self.save_sel_markers = False
            self.checked = True

        else:
            self.checked = True

        self.reset_gui(reloaded=True)

    def parse_markers(self):
        if self.save_asc_markers:
            active_event_marker, asc_phase_markers = \
                self.get_active_event_and_phase_markers()
            self.phase_markers = asc_phase_markers[:]
            asc_phase_markers.insert(0, active_event_marker)
            self.selected_markers = asc_phase_markers
            self.active_event = active_event_marker.get_event()
            # self.ae_name = self.active_event.name
            start_time = datetime.datetime(1970, 1, 1)
            event_time = self.active_event.time
            event_date = start_time.utcfromtimestamp(event_time)
            self.ae_name = event_date.strftime("%Y%m%d_%H%M%S%f")

        elif self.save_sel_markers:
            self.selected_markers = self.get_selected_markers()

        else:
            self.selected_markers = self.get_markers()

    def save_markers(self):
        if self.name_to_path and self.name_as_filename:
            save_path = self.root_folder + "/" + self.ae_name
            self.folder_check(save_path)
            self.save_path = save_path + "/" + self.ae_name

        elif self.name_to_path:
            save_path = self.root_folder + "/" + self.ae_name
            self.folder_check(save_path)
            self.save_path = save_path + "/picks"

        elif self.name_as_filename:
            self.save_path = self.root_folder + "/" + self.ae_name + "_picks"

        else:
            self.save_path = self.root_folder + "/picks"

        Marker.save_markers(self.selected_markers, self.save_path + '.txt')

    def save_velest(self):
        start_time = datetime.datetime(1970, 1, 1)
        event_time = self.active_event.time
        event_date = start_time.utcfromtimestamp(event_time)
        event_year = int(str(event_date.year)[2:])
        event_month = event_date.month
        event_day = event_date.day
        event_hour = event_date.hour
        event_min = event_date.minute
        event_sec = event_date.second + (
            round(event_date.microsecond / 10 ** 6, 2))
        event_lat = round(self.active_event.lat, 4)
        event_lon = round(self.active_event.lon, 4)
        event_depth = round(self.active_event.depth / 10 ** 3, 2)
        event_mag = 0.0  # round(self.active_event.magnitude, 2)

        if event_lat >= 0:
            event_NS = "N"
        else:
            event_NS = "S"

        if event_lon >= 0:
            event_EW = "E"
        else:
            event_EW = "W"

        event_line = "{:02d}{:02d}{:02d} {:02d}{:02d} {:05.2f} {:7.4f}{} \
{:8.4f}{}{:8.2f}{:7.2f}     99  0.0 0.00  1.0  1.0 \n".format(
            event_year, event_month, event_day, event_hour, event_min,
            event_sec, event_lat, event_NS, event_lon, event_EW, event_depth,
            event_mag
        )

        self.save_path_velest = self.save_path + "_velest.txt"

        if self.name_to_path or self.name_as_filename:
            velest_file = open(self.save_path_velest, 'w')
        else:
            velest_file = open(self.save_path_velest, 'a+')

        velest_file.write(event_line)

        phase_i = 0
        phase_num = len(self.phase_markers)

        for phase_marker in self.phase_markers:
            phase_i += 1
            phase_station = list(phase_marker.nslc_ids)[0][1]
            phase_name = phase_marker._phasename
            phase_tmin = phase_marker.tmin
            phase_tmax = phase_marker.tmax
            # phase_unc = phase_marker._uncertainty
            if phase_tmin == phase_tmax:
                phase_unc = 0.1
                phase_time = phase_tmin
            else:
                phase_unc = phase_tmax - phase_tmin
                phase_mid_int = (phase_tmax - phase_tmin) * 0.5
                phase_time = phase_tmin + phase_mid_int

            phase_rtime = round(phase_time - event_time, 2)

            if phase_unc <= 0.1:
                phase_unc_class = 0
            elif 0.1 < phase_unc <= 0.2:
                phase_unc_class = 1
            elif 0.2 < phase_unc <= 0.5:
                phase_unc_class = 2
            else:
                phase_unc_class = 3

            if phase_i % 6.0 == 0 or phase_i == phase_num:
                phase_block = "{:4s}{}{}{:06.2f}\n".format(phase_station,
                                                           phase_name, phase_unc_class, phase_rtime)
            else:
                phase_block = "{:4s}{}{}{:06.2f}".format(phase_station,
                                                         phase_name, phase_unc_class, phase_rtime)

            velest_file.write(phase_block)

        velest_file.write('\n')
        velest_file.close()

    def save_quakeml(self):
        from obspy.core.event import Catalog, Event, Pick, CreationInfo, \
            WaveformStreamID
        from obspy.core.event.base import QuantityError, Comment
        from obspy.core.event.origin import Origin
        from obspy import UTCDateTime

        creation_string = CreationInfo(creation_time=UTCDateTime.now())

        start_time = datetime.datetime(1970, 1, 1)
        event_time = self.active_event.time
        #event_date = UTCDateTime(start_time.utcfromtimestamp(event_time))
        #event_lat = round(self.active_event.lat, 4)
        #event_lon = round(self.active_event.lon, 4)
        #event_depth = round(self.active_event.depth / 10 ** 3, 2)
        #event_mag = 0.0  # round(self.active_event.magnitude, 2)
        latlon_error = QuantityError(uncertainty=0.009)  # equals 1 km

        event = Event()
        event.creation_info = creation_string

        for phase_marker in self.phase_markers:

            phase_station = list(phase_marker.nslc_ids)[0][1]
            phase_network = list(phase_marker.nslc_ids)[0][0]
            phase_channel = list(phase_marker.nslc_ids)[0][3]
            wav_id = WaveformStreamID(station_code=phase_station,
                                      channel_code=phase_channel,
                                      network_code=phase_network)
            phase_name = phase_marker.get_phasename()
            phase_tmin = phase_marker.tmin
            phase_tmax = phase_marker.tmax
            if phase_tmin == phase_tmax:
                phase_unc = 0.1
                phase_time = phase_tmin
            else:
                phase_unc = (phase_tmax - phase_tmin) * 0.5
                phase_mid_int = (phase_tmax - phase_tmin) * 0.5
                phase_time = phase_tmin + phase_mid_int

            phase_time_utc = UTCDateTime(start_time.utcfromtimestamp(phase_time))
            if phase_marker.get_polarity() == 1:
                phase_polarity = "positive"
            elif phase_marker.get_polarity() == -1:
                phase_polarity = "negative"
            else:
                phase_polarity = "undecidable"

            pick = Pick(time=phase_time_utc, phase_hint=phase_name, waveform_id=wav_id, polarity=phase_polarity,
                        time_errors=QuantityError(uncertainty=phase_unc), method_id="snuffler")
            event.picks.append(pick)

        # Add origin.
        # Use loc of station with earliest pick and a delay of -0.12 s, a ray travelling up from below at 250.0 m depth.
        earliest = min([p.time for p in event.picks])
        approx_event_time = earliest - 0.12
        stations = self.get_stations()
        if phase_station not in list([s.station for s in stations]):
            print("ERROR: Station %s not found in station file!" % phase_station)
        for station in stations:
            if station.station == phase_station:
                latitude = station.lat
                longitude = station.lon
        origin = Origin(time=approx_event_time, time_errors=QuantityError(uncertainty=0.12),
                        longitude=longitude, latitude_errors=latlon_error,
                        latitude=latitude, longitude_errors=latlon_error,
                        depth=250.0, depth_errors=QuantityError(uncertainty=250.0),
                        method_id="dummy")
        event.origins.append(origin)
        event.preferred_origin_id = str(event.origins[0].resource_id)

        cat = Catalog()
        cat.append(event)
        print(event)
        self.save_path_quakeml = self.save_path + ".xml"
        print(self.save_path_quakeml)
        cat.write(str(self.save_path_quakeml), format="QUAKEML")
        self.save_path_nlloc = self.save_path + ".obs"
        cat.write(str(self.save_path_nlloc), format="NLLOC_OBS")

    def call(self):
        self.check_params()

        if self.checked:
            self.parse_markers()
            self.save_markers()
            print('Markers saved to: {}.'.format(self.save_path + '.txt'))
            self.show_message('Success', ' Markers saved to: {}.'.format(
                self.save_path + '.txt'
            ))

            if self.velest_format:
                self.save_velest()
                print('Velest format saved to: {}.'.format(
                    self.save_path_velest
                ))
                self.show_message(
                    'Success', ' Velest format saved to: {}.'.format(
                        self.save_path_velest
                    ))
            if self.quakeml_format:
                self.save_quakeml()
                print('QuakeML format saved to: {}.'.format(
                    self.save_path_quakeml
                ))
                self.show_message(
                    'Success', ' QuakeML format saved to: {0} and nonlinloc saved to: {1}.'.format(
                        self.save_path_quakeml, self.save_path_nlloc
                    ))
            else:
                pass

        else:
            self.save_sel_markers = False
            self.save_asc_markers = False
            self.name_to_path = False
            self.name_as_filename = False
            self.velest_format = False
            self.quakeml_format = False
            self.reset_gui(reloaded=True)
            self.warn(
                ' Please select appropriate parameters. All markers are \
saved in case none is selected.\n\n\t\tSelection cleared.'
            )


def __snufflings__():
    """Returns a list of snufflings to be exported by this module."""

    return [QuickSave()]
