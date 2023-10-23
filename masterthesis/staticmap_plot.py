#!/usr/bin/env python3

import extract_post_gauges
import matplotlib.colors as colors
import matplotlib.figure as figure
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import numpy.typing as typing

CATCHMENT_FILE = (
    "/data/home/richi/master_thesis/model_MAR/staticmaps/wflow_subcatch.map"
)
DEM_FILE = "/data/home/richi/master_thesis/model_MAR/staticmaps/wflow_dem.map"
RIVER_FILE = "/data/home/richi/master_thesis/model_MAR/staticmaps/wflow_river.map"
GAUGE_FILE = "/data/home/richi/master_thesis/model_MAR/staticmaps/wflow_gauges.map"
STREAMORDER_FILE = (
    "/data/home/richi/master_thesis/model_MAR/staticmaps/wflow_streamorder.map"
)
LDD_FILE = "/data/home/richi/master_thesis/model_MAR/staticmaps/wflow_ldd.map"
LANDUSE_FILE = "/data/home/richi/master_thesis/model_MAR/staticmaps/wflow_landuse.map"

DEM_PGF = "/data/home/richi/master_thesis/thesis/images/wflow_dem.pgf"
RIVER_PGF = "/data/home/richi/master_thesis/thesis/images/wflow_river.pgf"
STREAMORDER_PGF = "/data/home/richi/master_thesis/thesis/images/wflow_streamorder.pgf"
LDD_PGF = "/data/home/richi/master_thesis/thesis/images/wflow_ldd.pgf"
LANDUSE_PGF = "/data/home/richi/master_thesis/thesis/images/wflow_landuse.pgf"
LANDUSE_LEVEL_1_PGF = (
    "/data/home/richi/master_thesis/thesis/images/wflow_landuse_level_1.pgf"
)

WIDTH_IN = extract_post_gauges.TEXTWITH_IN
HEIGHT_IN = WIDTH_IN

LANDUSE_CODES = {
    1: {
        "name": "Artificial Surfaces",
        1: {
            "name": "Urban fabric",
            1: {"name": "Continuous urban fabric"},
            2: {"name": "Discontinuous urban fabric"},
        },
        2: {
            "name": "Industrial, commercial and transport units",
            1: {"name": "Industrial or commercial units"},
            2: {"name": "Road and rail networks and associated land"},
            3: {"name": "Port areas"},
            4: {"name": "Airports"},
        },
        3: {
            "name": "Mine, dump and construction sites",
            1: {"name": "Mineral extraction sites"},
            2: {"name": "Dump sites"},
            3: {"name": "Construction sites"},
        },
        4: {
            "name": "Artificial, non-agricultural vegetated areas",
            1: {"name": "Green urban areas"},
            2: {"name": "Sport and leisure facilities"},
        },
    },
    2: {
        "name": "Agricultural Areas",
        1: {
            "name": "Arable land",
            1: {"name": "Non-irrigated arable land"},
            2: {"name": "Permanently irrigated land"},
            3: {"name": "Rice fields"},
        },
        2: {
            "name": "Permanent crops",
            1: {"name": "Vineyards"},
            2: {"name": "Fruit trees and berry plantations"},
            3: {"name": "Olive groves"},
        },
        3: {
            "name": "Pastures",
            1: {"name": "Pastures"},
        },
        4: {
            "name": " Heterogeneous agricultural areas",
            1: {"name": "Annual crops associated with permanent crops"},
            2: {"name": "Complex cultivation patterns"},
            3: {
                "name": "Land principally occupied by agriculture, with significant areas of natural vegetation"
            },
            4: {"name": "Agro-forestry areas"},
        },
    },
    3: {
        "name": "Forest and Seminatural Areas",
        1: {
            "name": "Forests",
            1: {"name": "Broad-leaved forest"},
            2: {"name": "Coniferous forest"},
            3: {"name": "Mixed forest"},
        },
        2: {
            "name": "Scrub and/or herbaceous associations",
            1: {"name": "Natural grassland"},
            2: {"name": "Moors and heathland"},
            3: {"name": "Sclerophyllous vegetation"},
            4: {"name": "Transitional woodland-scrub"},
        },
        3: {
            "name": "Open spaces with little or no vegetation",
            1: {"name": "Beaches, dunes, sands"},
            2: {"name": "Bare rocks"},
            3: {"name": "Sparsely vegetated areas"},
            4: {"name": "Burnt areas"},
            5: {"name": "Glaciers and perpetual snow"},
        },
    },
    4: {
        "name": "Wetlands",
        1: {
            "name": "Inland wetlands",
            1: {"name": "Inland marshes"},
            2: {"name": "Peat bogs"},
        },
        2: {
            "name": "Marine wetlands",
            1: {"name": "Salt marshes"},
            2: {"name": "Salines"},
            3: {"name": "Intertidal flats"},
        },
    },
    5: {
        "name": "Water Bodies",
        1: {
            "name": "Inland waters",
            1: {"name": "Water courses"},
            2: {"name": "Water bodies"},
        },
        2: {
            "name": "Marine waters",
            1: {"name": "Coastal lagoons"},
            2: {"name": "Estuaries"},
            3: {"name": "Sea and ocean"},
        },
    },
}


def savefig(fig: figure.Figure, filename: str) -> None:
    fig.tight_layout()
    if filename.endswith(".pgf"):
        fig.savefig(filename, backend="pgf")
    else:
        fig.savefig(filename)


def getCatchmentMask(catchment_file: str) -> typing.NDArray:
    catchment_mask = extract_post_gauges.getRaster(catchment_file)
    catchment_mask[catchment_mask == 1] = 0
    catchment_mask[catchment_mask == np.nan] = 1
    return catchment_mask


def generateDem(dem_file: str, catchment_file: str) -> typing.NDArray:
    return ma.masked_array(
        extract_post_gauges.getRaster(dem_file), mask=getCatchmentMask(catchment_file)
    )


def writeDem(dem_file: str, catchment_file: str, outfile: str) -> None:
    fig, ax = plt.subplots()
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    dem = ma.masked_array(
        extract_post_gauges.getRaster(dem_file), mask=getCatchmentMask(catchment_file)
    )
    dem_ax = ax.imshow(dem)
    label = "Elevation [\\si{\\meter}]" if outfile.endswith(".pgf") else "Elevation [m]"
    fig.colorbar(
        dem_ax,
        ax=ax,
        label=label,
        orientation="horizontal",
    )

    savefig(fig, outfile)
    plt.close(fig)


def writeRiverAndGauges(
    dem_file: str, catchment_file: str, river_file: str, gauges_file: str, outfile: str
) -> None:
    river = extract_post_gauges.getRaster(river_file, no_data=-1)
    river_mask = np.zeros(river.shape)
    river_mask[river == -1] = 1
    river_masked = ma.masked_array(river, mask=river_mask)

    gauges = extract_post_gauges.getRaster(gauges_file, no_data=-1)
    gauge_mask = np.zeros(gauges.shape)
    gauge_mask[gauges == -1] = 1
    gauges_masked = ma.masked_array(gauges, mask=gauge_mask)

    fig, ax = plt.subplots()
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    dem_ax = ax.imshow(generateDem(dem_file, catchment_file), alpha=0.6)
    label = "Elevation [\\si{\\meter}]" if outfile.endswith(".pgf") else "Elevation [m]"
    fig.colorbar(dem_ax, ax=ax, label=label, orientation="horizontal")
    river_map = colors.ListedColormap(["blue"])
    ax.imshow(river_masked, cmap=river_map, interpolation="none")
    gauge_colors = ["white" for _ in range(9)]
    gauge_colors[2] = "red"
    gauge_map = colors.ListedColormap(gauge_colors)
    ax.imshow(gauges_masked, cmap=gauge_map, interpolation="none")

    savefig(fig, outfile)
    plt.close(fig)


def writeStreamOrder(catchment_file: str, streamorder_file: str, outfile: str) -> None:
    streamorder = ma.masked_array(
        extract_post_gauges.getRaster(streamorder_file),
        mask=getCatchmentMask(catchment_file),
    )

    fig, ax = plt.subplots()
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    stream_ax = ax.imshow(streamorder, cmap="plasma", interpolation="none")
    cbar = fig.colorbar(stream_ax, ax=ax, label="Streamorder", orientation="horizontal")
    ticks = range(int(np.nanmin(streamorder)), int(np.nanmax(streamorder)) + 1)
    cbar.set_ticks(
        ticks=ticks,
        labels=[str(i) for i in ticks],
    )

    savefig(fig, outfile)
    plt.close(fig)


def writeLdd(catchment_file: str, ldd_file: str, outfile: str) -> None:
    ldd = ma.masked_array(
        extract_post_gauges.getRaster(ldd_file, -1),
        mask=getCatchmentMask(catchment_file),
    )

    fig, ax = plt.subplots()
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    ldd_ax = ax.imshow(ldd, cmap="viridis", interpolation="none")
    cbar = fig.colorbar(ldd_ax, ax=ax, label="Direction", orientation="horizontal")
    ticks = range(int(np.nanmin(ldd)), int(np.nanmax(ldd)) + 1)
    cbar.set_ticks(
        ticks=ticks,
        labels=[str(i) for i in ticks],
    )

    savefig(fig, outfile)
    plt.close(fig)


def writeLandUse(
    catchment_file: str, landuse_file: str, outfile: str, outfile_level_1: str = ""
) -> None:
    landuse = ma.masked_array(
        extract_post_gauges.getRaster(landuse_file, -1),
        mask=getCatchmentMask(catchment_file),
    )

    fig, ax = plt.subplots()
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    landuse_ax = ax.imshow(landuse, cmap="viridis", interpolation="none")
    fig.colorbar(landuse_ax, ax=ax, label="Landuse ID", orientation="horizontal")

    savefig(fig, outfile)
    plt.close(fig)

    if outfile_level_1:
        level1 = landuse // 100
        level1_colors = {
            1: "black",
            2: "yellow",
            3: "green",
            4: "mediumaquamarine",
            5: "blue",
        }
        used_ids = np.unique(level1.compressed())
        level1_map = colors.ListedColormap([level1_colors[i] for i in used_ids])
        patches = [
            mpatches.Patch(color=level1_colors[i], label=LANDUSE_CODES[i]["name"])
            for i in used_ids
        ]

        fig, ax = plt.subplots()
        fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN)
        ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

        ax.imshow(level1, cmap=level1_map, interpolation="none")
        fig.legend(handles=patches, loc="lower left")

        savefig(fig, outfile_level_1)
        plt.close(fig)


def main() -> None:
    plt.rcParams.update(
        {
            "pgf.preamble": r"\usepackage{siunitx}",
            "font.family": "serif",
            "pgf.rcfonts": False,
        }
    )

    writeDem(DEM_FILE, CATCHMENT_FILE, DEM_PGF)
    writeRiverAndGauges(DEM_FILE, CATCHMENT_FILE, RIVER_FILE, GAUGE_FILE, RIVER_PGF)
    writeStreamOrder(CATCHMENT_FILE, STREAMORDER_FILE, STREAMORDER_PGF)
    writeLdd(CATCHMENT_FILE, LDD_FILE, LDD_PGF)
    writeLandUse(CATCHMENT_FILE, LANDUSE_FILE, LANDUSE_PGF, LANDUSE_LEVEL_1_PGF)


if __name__ == "__main__":
    main()
