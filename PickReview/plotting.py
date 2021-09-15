import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from obspy.core.util.attribdict import AttribDict
from obspy.io.nlloc.util import read_nlloc_scatter
from obspy import read_events
import numpy as np
from pyproj import Proj
import pandas as pd
from obspy.geodetics.base import degrees2kilometers
from obspy import Catalog

plt.rcParams.update({'font.size': 14})

# Get station information
nllstafile = "/home/genevieve/research/PickReview/NonLinLoc/time_tables/gtsrc_stations_welev.dat"
stalst = pd.read_csv(nllstafile, delimiter=" ", skiprows=0,
                     usecols=[1, 3, 4, 5, 6],
                     header=None,
                     names=["station", "latitude", "longitude", "depth", "elevation"])

is_surface = stalst["depth"] == 0
stalst = stalst[is_surface]
stalats = stalst["latitude"].to_list()
stalons = stalst["longitude"].to_list()
staname = stalst["station"].to_list()

# CaMI.FRS Wells
injwell = AttribDict({"name": "injection", "lat": 50.45041, "lon": -112.12067, "color": "m"})
geowell = AttribDict({"name": "geophysics", "lat": 50.45031, "lon": -112.12087, "color": "b"})
chemwell = AttribDict({"name": "geochemistry", "lat": 50.45061, "lon": -112.12036, "color": "g"})


def plot_scatter(hypfile, figfile=None, picked_stations=[]):
    cat = read_events(hypfile, format="NLLOC_HYP")

    # lat, long to X,Y
    p = Proj(proj='eqc', lat_0=injwell.lat, lon_0=injwell.lon)
    evlat = cat[0].preferred_origin().latitude
    evlon = cat[0].preferred_origin().longitude
    xev, yev = p(evlon, evlat)
    xev *= 1e-3
    yev *= 1e-3
    zev = cat[0].preferred_origin().depth * 1e-3
    print(f"lat = {evlat}, long={evlon}, x = {xev}, y={yev}, z={zev}")

    scatter_cloud = read_nlloc_scatter(hypfile.replace(".hyp", ".scat"))
    scatter_cloud = sorted(scatter_cloud, key=lambda x: x[3])

    x = [s[0] for s in scatter_cloud]
    y = [s[1] for s in scatter_cloud]
    z = [s[2] for s in scatter_cloud]
    pdf = [s[3] for s in scatter_cloud]

    xbest = x[-1]
    ybest = y[-1]
    zbest = z[-1]

    stalst = pd.read_csv("/home/genevieve.savard/cami/catalog/stations_XYZ.csv")
    xs = stalst["x"].values
    ys = stalst["y"].values
    zs = stalst["z"].values
    xsbh = np.append(stalst.loc[stalst.station.str.startswith("BH")].x.values, 0)
    ysbh = np.append(stalst.loc[stalst.station.str.startswith("BH")].y.values, 0)
    zsbh = np.append(stalst.loc[stalst.station.str.startswith("BH")].z.values, 0)
    xspick = []
    yspick = []
    zspick = []
    for sta in picked_stations:
        if "FRS" in sta:
            continue
        xspick.append(stalst.loc[stalst.station == sta].x.values[0])
        yspick.append(stalst.loc[stalst.station == sta].y.values[0])
        zspick.append(stalst.loc[stalst.station == sta].z.values[0])

    fig, axs = plt.subplots(2, 2, figsize=(14, 14))
    axs[1][1].set_axis_off()

    # X-Y
    axs[0][0].scatter(x, y, c=pdf)
    axs[0][0].scatter(xbest, ybest, s=150, c="r", marker="*", edgecolors="k", linewidths=1.0)
    # axs[0][0].scatter(xev, yev, s=50, c="c", marker="*")
    axs[0][0].scatter(xs, ys, s=50, c="c", marker="s", edgecolors="k", linewidths=1.0, zorder=40)
    axs[0][0].scatter(xspick, yspick, s=50, c="m", marker="s", edgecolors="k", linewidths=1.0, zorder=41)
    axs[0][0].plot(0, 0, marker="o", markersize=15, color="m", markeredgecolor="k", markeredgewidth=4)  # Injection well
    axs[0][0].set_xlabel("Easting (km)")
    axs[0][0].set_ylabel("Northing (km)")
    axs[0][0].xaxis.set_label_position('top')

    # X-Z
    # geology
    axs[1][0].axhline(0, -1, 1, color="k", zorder=1)
    axs[1][0].axhspan(0.025, 0.099, facecolor="brown", alpha=0.25, zorder=1, ec="k", lw=1)  # bedrock/Dinosaur Park top
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
    axs[1][0].scatter(xspick, zspick, s=50, c="m", marker="s", edgecolors="k", linewidths=1.0, zorder=41)
    axs[1][0].plot(0, 0, marker="o", markersize=15, color="m", markeredgecolor="k", markeredgewidth=4)  # Injection well
    axs[1][0].invert_yaxis()
    axs[1][0].set_ylabel("Depth (km)")
    # ax.set_yticklabels([])
    # ax.set_xticklabels([])

    # Y-Z
    # geology
    axs[0][1].axvline(0, -1, 1, color="k", zorder=1)
    axs[0][1].axvspan(0.025, 0.099, facecolor="brown", alpha=0.25, zorder=1, ec="k", lw=1)  # bedrock/Dinosaur Park top
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
    axs[0][1].scatter(zspick, yspick, s=50, c="m", marker="s", edgecolors="k", linewidths=1.0, zorder=41)
    axs[0][1].plot(0, 0, marker="o", markersize=15, color="m", markeredgecolor="k", markeredgewidth=4)  # Injection well
    axs[0][1].set_yticklabels([])
    axs[0][1].set_xlabel("Depth (km)")

    plt.subplots_adjust(wspace=0.01, hspace=0.01)
    plt.show()
    if figfile:
        plt.savefig(figfile)

    plt.close()


