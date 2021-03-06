import datetime
import obspy.core.event
from pyrocko.util import str_to_time
from pyrocko.gui.marker import PhaseMarker, Marker
from obspy.core.event import Catalog, Event, Pick, CreationInfo, \
    WaveformStreamID
from obspy.core.event.base import QuantityError, Comment
from obspy.core.event.origin import Origin
from obspy import UTCDateTime
from pyrocko import obspy_compat
import logging

obspy_compat.plant()
Logger = logging.getLogger(__name__)

# Parameters to define a dummy origin when creating Obspy Event object from picked phases
DEFAULT_LATLONG_ERROR = QuantityError(uncertainty=0.009)  # Default Lat/long origin error for event creation: equals 1 km
DEFAULT_TIME_ERROR = QuantityError(uncertainty=0.12)  # Default origin time error
REFERENCE_UTC_TIME = datetime.datetime(1970, 1, 1)
DEFAULT_DEPTH_ERROR = QuantityError(uncertainty=250.0)
DEFAULT_DEPTH = 250.0  # Default Origin depth
DEFAULT_LATITUDE = 50.45041  # Use injection well
DEFAULT_LONGITUDE = -112.12067  # Use injection well
DEFAULT_VP = 2000.
DEFAULT_TRAVELTIME_BELOW_STATION = DEFAULT_DEPTH / DEFAULT_VP
DEFAULT_PICK_UNCERTAINTY = 0.06


def picks2markers(picks, event=None, phase=True, kinds=(1, 2)):
    """
    Convert list of Obspy Picks to Pyrocko Markers list
    :param event: Obspy Event
    :param phase: Make Phase Markers
    :param kinds: tuple of length 2 with Marker "kind" (color) for P, then S
    :param picks: list of Obspy Pick
    :return: list of Pyrocko PhaseMarker
    """
    if phase and not event:
        Logger.error("Must supply parameter event if phase=True")
    # Convert event Obspy -> Pyrocko
    cat = obspy.Catalog(events=[event])
    evpyro = obspy_compat.to_pyrocko_events(cat)[0]

    # Make marker list
    markers = []
    for p in picks:
        if not p.waveform_id.location_code:
            nscl = (p.waveform_id.network_code, p.waveform_id.station_code, "", p.waveform_id.channel_code)
        else:
            nscl = (p.waveform_id.network_code, p.waveform_id.station_code, p.waveform_id.location_code, p.waveform_id.channel_code)
        kind = kinds[0] if "p" in p.phase_hint.lower() else kinds[1]
        pick_time = str_to_time(p.time.strftime("%Y-%m-%d %H:%M:") + "%f" % (p.time.second + p.time.microsecond * 1e-6))

        if phase:
            m = PhaseMarker(nslc_ids=[nscl], tmin=pick_time, tmax=pick_time, kind=kind, event=evpyro, phasename=p.phase_hint)
        else:
            m = Marker(nslc_ids=[nscl], tmin=pick_time, tmax=pick_time, kind=kind)
        markers.append(m)
    return markers


def fix_picks_ids(event, stream, method=None):
    """
    Fix the picks in an obspy Event don't have the matching waveform ID with the corresponding Obspy stream
     Add network and dummy channel if station info available. N or 1 used as default for S phases, Z for P phases
    :param method: String, only include pick for which string is in Pick.method_id.id
    :param event: Obspy.Event
    :param stream: Obspy.Stream
    :return: list of Obspy Pick, Obspy.Event
    """
    newpicks = []  # Will create new pick list because in-place modification don't always work...
    for p in event.picks:
        newp = p.copy()
        if not p.waveform_id.network_code or not p.waveform_id.channel_code:
            station = p.waveform_id.station_code
            traces = stream.select(station=station)
            if not traces:
                Logger.warning(f"No waveforms in stream matching station {station} found. Remove those picks.")
                continue
            network = traces[0].stats.network
            if method:
                if method.lower() not in p.method_id.id.lower():
                    Logger.info(f"Skipping pick with method_id = {p.method_id}.")
                    continue
            if "s" in p.phase_hint.lower():  # Choose N or 1 as channel
                channel = traces.select(channel="*[N1]")[0].stats.channel
            elif "p" in p.phase_hint.lower():  # Choose Z channel
                channel = traces.select(channel="*Z")[0].stats.channel
            newp.waveform_id = WaveformStreamID(station_code=station, network_code=network, channel_code=channel)
            newpicks.append(newp)
        else:
            newpicks.append(newp)

    return newpicks, obspy.core.event.Event(picks=newpicks)


def markers_to_quakeml(phase_markers, picks_only=True, inv=None, origin=None):
    """
    Convert list of Pyrocko Phase Markers to Obspy event
    :param phase_markers: list of Pyrocko Phase Markers
    :param picks_only: If creating Obspy Event with only picks, no origin
    :param inv: Obspy Inventory with station info
    :param origin: Obspy Origin if using origin info from previous location run
    :return: Obspy Event
    """
    event = Event()

    # Creation info
    creation_string = CreationInfo(creation_time=UTCDateTime.now())
    event.creation_info = creation_string

    # Loop over phase markers and convert to Obspy Pick
    stations = []
    sta_time_tup = []
    for phase_marker in phase_markers:
        if not isinstance(phase_marker, PhaseMarker):
            Logger.info("Skipping marker that's not a PhaseMarker")
            continue
        phase_station = list(phase_marker.nslc_ids)[0][1]
        stations.append(phase_station)
        phase_network = list(phase_marker.nslc_ids)[0][0]
        phase_channel = list(phase_marker.nslc_ids)[0][3]
        wav_id = WaveformStreamID(station_code=phase_station,
                                  channel_code=phase_channel,
                                  location_code="",
                                  network_code=phase_network)
        phase_name = phase_marker.get_phasename()  # _phasename

        phase_tmin = phase_marker.tmin
        phase_tmax = phase_marker.tmax
        if phase_tmin == phase_tmax:
            phase_unc = DEFAULT_PICK_UNCERTAINTY
            phase_time = phase_tmin
        else:
            phase_unc = (phase_tmax - phase_tmin) * 0.5
            phase_mid_int = (phase_tmax - phase_tmin) * 0.5
            phase_time = phase_tmin + phase_mid_int

        phase_time_utc = UTCDateTime(REFERENCE_UTC_TIME.utcfromtimestamp(phase_time))

        if phase_marker.get_polarity() == 1:
            phase_polarity = "positive"
        elif phase_marker.get_polarity() == -1:
            phase_polarity = "negative"
        else:
            phase_polarity = "undecidable"

        pick = Pick(time=phase_time_utc, phase_hint=phase_name, waveform_id=wav_id, polarity=phase_polarity,
                    time_errors=QuantityError(uncertainty=phase_unc), method_id="snuffler")
        event.picks.append(pick)
        sta_time_tup.append((wav_id.get_seed_string(), pick.time))

    if not picks_only:
        # Add origin.

        if not origin:
            # Use loc of station with earliest pick and dummy origin situated below it at DEFAULT_DEPTH
            sta_time_tup.sort(key=lambda x: x[1], reverse=False)
            earliest_sta_id = sta_time_tup[0][0]
            dummy_origin_time = sta_time_tup[0][1] - DEFAULT_TRAVELTIME_BELOW_STATION

            if inv:
                try:
                    latitude = inv.get_coordinates(earliest_sta_id)["latitude"]
                    longitude = inv.get_coordinates(earliest_sta_id)["longitude"]
                except:
                    latitude = DEFAULT_LATITUDE
                    longitude = DEFAULT_LONGITUDE
            else:
                latitude = DEFAULT_LATITUDE
                longitude = DEFAULT_LONGITUDE

            origin = Origin(time=dummy_origin_time, time_errors=DEFAULT_TIME_ERROR,
                            longitude=longitude, latitude_errors=DEFAULT_LATLONG_ERROR,
                            latitude=latitude, longitude_errors=DEFAULT_LATLONG_ERROR,
                            depth=DEFAULT_DEPTH, depth_errors=DEFAULT_DEPTH_ERROR,
                            method_id="dummy")
        event.origins.append(origin)
        event.preferred_origin_id = str(event.origins[0].resource_id)

    return event


def write_obs(phase_markers, fname):
    """
    :param phase_markers: list of Pyrocko Phase Markers
    :param fname: output filename for NLLOC OBS phase file
    """
    with open(fname, "w") as outfile:
        for phase_marker in phase_markers:
            if not isinstance(phase_marker, PhaseMarker):
                Logger.info("Skipping marker that isn't a PhaseMarker")
                continue
            station = list(phase_marker.nslc_ids)[0][1]
            channel = list(phase_marker.nslc_ids)[0][3]
            phase_name = phase_marker.get_phasename()  # _phasename

            # Get time and uncertainty
            phase_tmin = phase_marker.tmin
            phase_tmax = phase_marker.tmax
            if phase_tmin == phase_tmax:
                phase_unc = DEFAULT_PICK_UNCERTAINTY
                phase_time = phase_tmin
            else:
                phase_unc = (phase_tmax - phase_tmin) * 0.5
                phase_mid_int = (phase_tmax - phase_tmin) * 0.5
                phase_time = phase_tmin + phase_mid_int
            phase_time_utc = UTCDateTime(REFERENCE_UTC_TIME.utcfromtimestamp(phase_time))

            # First motion polarity
            if phase_marker.get_polarity() == 1:
                polarity = "+"
            elif phase_marker.get_polarity() == -1:
                polarity = "-"
            else:
                polarity = "?"

            # Other defaults
            instrument = "?"
            onset = "?"
            date = phase_time_utc.strftime("%Y%m%d")
            hourmin = phase_time_utc.strftime("%H%M")
            seconds = phase_time_utc.strftime("%S.%f")[:-2]
            error_type = "GAU"
            coda_duration = -1
            amplitude = -1  # Maximum peak-to-peak amplitude
            period = -1  # Period of amplitude reading
            priorwt = 1

            # Write phase to file
            line = f"{station:6s} " \
                   f"{instrument:4s} " \
                   f"{channel:4s} " \
                   f"{onset:1s} " \
                   f"{phase_name:6s} " \
                   f"{polarity:1s} " \
                   f"{date:6s} " \
                   f"{hourmin:4s} " \
                   f"{seconds:7s} " \
                   f"{error_type:4s} " \
                   f"{phase_unc:9.2e} " \
                   f"{coda_duration:9.2e} " \
                   f"{amplitude:9.2e} " \
                   f"{period:9.2e} " \
                   f"{priorwt:9.2e}\n"
            outfile.write(line)
