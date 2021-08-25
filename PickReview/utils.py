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
obspy_compat.plant()

import logging
Logger = logging.getLogger(__name__)


def markers_to_quakeml(phase_markers, inv):
    """
    Convert list of Pyrocko Phase Markers to Obspy event
    :param phase_markers: list of Pyrocko Phase Markers
    :param inv: Obspy Inventory with station info
    :return: Obspy Event
    """
    creation_string = CreationInfo(creation_time=UTCDateTime.now())

    start_time = datetime.datetime(1970, 1, 1)
    latlon_error = QuantityError(uncertainty=0.009)  # equals 1 km

    event = Event()
    event.creation_info = creation_string

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
        phase_name = phase_marker.get_phasename() #_phasename

        phase_tmin = phase_marker.tmin
        phase_tmax = phase_marker.tmax
        if phase_tmin == phase_tmax:
            phase_unc = 0.05
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
        sta_time_tup.append((wav_id.get_seed_string(), pick.time))

    # Add origin.
    # Use loc of station with earliest pick and a delay of -0.12 s, a ray travelling up from below at 250.0 m depth.
    sta_time_tup.sort(key=lambda x: x[1], reverse=False)
    earliest_sta_id = sta_time_tup[0][0]
    approx_event_time = sta_time_tup[0][1] - 0.12

    try:
        latitude = inv.get_coordinates(earliest_sta_id)["latitude"]
        longitude = inv.get_coordinates(earliest_sta_id)["longitude"]
    except:
        latitude = 50.45041  # Use injection well
        longitude = -112.12067

    origin = Origin(time=approx_event_time, time_errors=QuantityError(uncertainty=0.12),
                    longitude=longitude, latitude_errors=latlon_error,
                    latitude=latitude, longitude_errors=latlon_error,
                    depth=250.0, depth_errors=QuantityError(uncertainty=250.0),
                    method_id="dummy")
    event.origins.append(origin)
    event.preferred_origin_id = str(event.origins[0].resource_id)
    cat = Catalog()
    cat.append(event)

    return cat


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
    evpyro = cat.to_pyrocko_events(cat)[0]

    # Make marker list
    markers = []
    for p in picks:
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
     Add network and dummy channel if station info available
    :param method: String, only include pick for which string is in Pick.method_id.id
    :param event: Obspy.Event
    :param stream: Obspy.Stream
    :return: list of Obspy Pick, Obspy.Event
    """
    newpicks = []  # Will create new pick list because in-place modification don't always work...
    for p in event.picks:
        station = p.waveform_id.station_code
        traces = stream.select(station=station)
        if not traces:
            Logger.warning(f"No waveforms in stream matching station {station} found. Remove those picks.")
            continue
        network = traces[0].stats.network
        if method:
            if method.lower() not in p.method_id.id.lower():
                Logger.info(f"Skipping pick with method_id = {p.method_id}.")
        if "s" in p.phase_hint.lower():  # Choose N or 1 as channel
            channel = traces.select(channel="*[1N]")[0].stats.channel
        elif "p" in p.phase_hint.lower():
            channel = traces.select(channel="*Z")[0].stats.channel
        newp = p.copy()
        newp.waveform_id = WaveformStreamID(station_code=station, network_code=network, channel_code=channel)
        newpicks.append(newp)

    return newpicks, obspy.core.event.Event(picks=newpicks)
