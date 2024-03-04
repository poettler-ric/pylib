#!/usr/bin/env python3

import matplotlib.colors as colors
import matplotlib.figure as figure
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import numpy.typing as typing
from osgeo import gdal, osr

import extract_post_gauges

PIXEL_WIDTH = 400
PIXEL_HEIGHT = -400
ORIGIN_X = 497000.0
ORIGIN_Y = 5295000.0
CRS = "EPSG:32633"
CRS_EPSG = 32633


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

DEM_OUTFILE = "/data/home/richi/master_thesis/thesis/images/wflow_dem.pgf"
RIVER_OUTFILE = "/data/home/richi/master_thesis/thesis/images/wflow_river.pgf"
STREAMORDER_OUTFILE = (
    "/data/home/richi/master_thesis/thesis/images/wflow_streamorder.png"
)
LDD_OUTFILE = "/data/home/richi/master_thesis/thesis/images/wflow_ldd.png"
LANDUSE_OUTFILE = "/data/home/richi/master_thesis/thesis/images/wflow_landuse.png"
LANDUSE_LEVEL_1_OUTFILE = (
    "/data/home/richi/master_thesis/thesis/images/wflow_landuse_level_1.pgf"
)

WIDTH_IN = extract_post_gauges.TEXTWITH_IN
HEIGHT_IN = WIDTH_IN
LAYOUT = "tight"

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
                "name": "Land principally occupied by agriculture"
                # "name": "Land principally occupied by agriculture, with significant areas of natural vegetation"
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


def getNameForLanduse(code):
    level1 = code // 100
    level2 = (code // 10) % 10
    level3 = code % 10
    return LANDUSE_CODES[level1][level2][level3]["name"]


def savefig(fig: figure.Figure, filename: str) -> None:
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
    fig, ax = plt.subplots(layout=LAYOUT)
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    dem = ma.masked_array(
        extract_post_gauges.getRaster(dem_file), mask=getCatchmentMask(catchment_file)
    )
    dem_ax = ax.imshow(dem)
    label = "Elevation [m.a.s.l.]"
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

    fig, ax = plt.subplots(layout=LAYOUT)
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    dem_ax = ax.imshow(generateDem(dem_file, catchment_file), alpha=0.6)
    label = "Elevation [m.a.s.l.]"
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

    unique_ids = np.sort(np.unique(streamorder[~streamorder.mask]))
    palette = list(
        generateColorMap(
            len(unique_ids),
            0,
            0,
            0,
            True,
        )
    )
    palette.reverse()
    colormap = colors.ListedColormap(palette)
    color_patches = [
        mpatches.Patch(color=color, label=f"Streamorder {streamorder_id:g}")
        for streamorder_id, color in zip(unique_ids, palette)
    ]

    fig, ax = plt.subplots(layout=LAYOUT)
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN * 2)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    stream_ax = ax.imshow(streamorder, cmap=colormap, interpolation="none")
    fig.legend(handles=color_patches, loc="lower left")

    savefig(fig, outfile)
    plt.close(fig)


def writeLdd(catchment_file: str, ldd_file: str, outfile: str) -> None:
    ldd = ma.masked_array(
        extract_post_gauges.getRaster(ldd_file, -1),
        mask=getCatchmentMask(catchment_file),
    )

    ldd_legends = {
        1: "Drainage direction to the southwest",
        2: "Drainage direction to the south",
        3: "Drainage direction to the southeast",
        4: "Drainage direction to the west",
        5: "No drainage direction (e.g. a pit)",
        6: "Drainage direction to the east",
        7: "Drainage direction to the northwest",
        8: "Drainage direction to the north",
        9: "Drainage direction to the northeast",
    }
    palette = list(
        generateColorMap(
            9,
            0,
            0,
            0,
            True,
        )
    )
    colormap = colors.ListedColormap(palette)
    color_patches = [
        mpatches.Patch(color=color, label=ldd_legends[direction])
        for direction, color in zip(range(1, 10), palette)
    ]

    fig, ax = plt.subplots(layout=LAYOUT)
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN * 2)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    ldd_ax = ax.imshow(ldd, cmap=colormap, interpolation="none")
    fig.legend(handles=color_patches, loc="lower left")
    # cbar = fig.colorbar(ldd_ax, ax=ax, label="Direction", orientation="horizontal")
    # ticks = range(int(np.nanmin(ldd)), int(np.nanmax(ldd)) + 1)
    # cbar.set_ticks(
    #     ticks=ticks,
    #     labels=[str(i) for i in ticks],
    # )

    savefig(fig, outfile)
    plt.close(fig)


def generateColorMap(count, red, green, blue, exclude_white=False):
    temp_count = count
    if exclude_white:
        temp_count += 1

    values = np.linspace(255, 0, temp_count)
    if exclude_white:
        values = values[:-1]

    red_scale = np.linspace(red, 255, temp_count)
    green_scale = np.linspace(green, 255, temp_count)
    blue_scale = np.linspace(blue, 255, temp_count)

    if exclude_white:
        red_scale = red_scale[:-1]
        green_scale = green_scale[:-1]
        blue_scale = blue_scale[:-1]

    return (
        (r, g, b)
        for r, g, b in zip(red_scale, green_scale, blue_scale)
        # (red if red else i, green if green else i, blue if blue else i)
        # for i in values
        # (red, green, blue, i)
        # for i in values
    )


def generateRedBlueMap(count):
    red = np.linspace(0, 1, count)
    blue = np.flip(red)
    return ((r, 0, b) for r, b in zip(red, blue))


def writeLandUse(
    catchment_file: str, landuse_file: str, outfile: str, outfile_level_1: str = ""
) -> None:
    landuse = ma.masked_array(
        extract_post_gauges.getRaster(landuse_file, -1),
        mask=getCatchmentMask(catchment_file),
    )

    level1_colors = {
        # 1: "black",
        1: (0, 0, 0),
        # 2: "yellow",
        2: (1, 1, 0),
        # 3: "green",
        3: (0, 0.5, 0),
        # 4: "mediumaquamarine",
        # 5: "blue",
    }

    unique_ids = np.sort(np.unique(landuse[~landuse.mask]))
    cathegories = {}
    # cathegorize ids
    for i in unique_ids:
        cathegory = i // 100
        if cathegory not in cathegories:
            cathegories[cathegory] = []
        cathegories[cathegory].append(i)

    palette = []
    for cathegory, values in cathegories.items():
        palette += list(
            generateColorMap(
                len(values),
                level1_colors[cathegory][0],
                level1_colors[cathegory][1],
                level1_colors[cathegory][2],
                True,
            )
        )
    colormap = colors.ListedColormap(palette)
    color_patches = [
        mpatches.Patch(color=color, label=getNameForLanduse(code))
        for code, color in zip(unique_ids, palette)
    ]

    fig, ax = plt.subplots(layout=LAYOUT)
    fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN * 3)
    ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    # landuse_ax = ax.imshow(
    #     landuse, cmap=colormap, interpolation="none", extent=(-1.5, 1, 1, 0)
    # )
    landuse_ax = ax.imshow(landuse, cmap=colormap, interpolation="none")
    fig.legend(handles=color_patches, loc="lower left")

    savefig(fig, outfile)
    plt.close(fig)

    if outfile_level_1:
        level1 = landuse // 100
        used_ids = np.unique(level1.compressed())
        level1_map = colors.ListedColormap([level1_colors[i] for i in used_ids])
        level1_patches = [
            mpatches.Patch(color=level1_colors[i], label=LANDUSE_CODES[i]["name"])
            for i in used_ids
        ]

        fig, ax = plt.subplots(layout=LAYOUT)
        fig.set_size_inches(w=WIDTH_IN, h=HEIGHT_IN)
        ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

        ax.imshow(level1, cmap=level1_map, interpolation="none")
        fig.legend(handles=level1_patches, loc="lower left")

        savefig(fig, outfile_level_1)
        plt.close(fig)


def create_raster(
    file_name,
    raster_array,
    origin,
    epsg=4326,
    pixel_width=10,
    pixel_height=10,
    nan_value=-9999.0,
    rdtype=gdal.GDT_Float32,
):
    """
    Convert a numpy.array to a GeoTIFF raster with the following parameters
    :param file_name: STR of target file name, including directory; must end on ".tif"
    :param raster_array: np.array of values to rasterize
    :param origin: TUPLE of (x, y) origin coordinates
    :param epsg: INT of EPSG:XXXX projection to use - default=4326
    :param pixel_height: INT of pixel height (multiple of unit defined with the EPSG number) - default=10m
    :param pixel_width: INT of pixel width (multiple of unit defined with the EPSG number) - default=10m
    :param nan_value: INT/FLOAT no-data value to be used in the raster (replaces non-numeric and np.nan in array)
                        default=-9999.0
    :param rdtype: gdal.GDALDataType raster data type - default=gdal.GDT_Float32 (32 bit floating point)
    """
    # check out driver
    driver = gdal.GetDriverByName("GTiff")

    # create raster dataset with number of cols and rows of the input array
    cols = raster_array.shape[1]
    rows = raster_array.shape[0]
    new_raster = driver.Create(file_name, cols, rows, 1, eType=rdtype)

    # apply geo-origin and pixel dimensions
    origin_x = origin[0]
    origin_y = origin[1]
    new_raster.SetGeoTransform((origin_x, pixel_width, 0, origin_y, 0, pixel_height))

    # replace np.nan values
    # raster_array[np.isnan(raster_array)] = nan_value

    # retrieve band number 1
    band = new_raster.GetRasterBand(1)
    band.SetNoDataValue(nan_value)
    band.WriteArray(raster_array)
    band.SetScale(1.0)

    # create projection and assign to raster
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    new_raster.SetProjection(srs.ExportToWkt())

    # release raster band
    band.FlushCache()


def prepareQgis():
    qgis_path = "/data/home/richi/master_thesis/thesis/qgis/catchment_maps"

    # dem
    dem = ma.masked_array(
        extract_post_gauges.getRaster(DEM_FILE), mask=getCatchmentMask(CATCHMENT_FILE)
    )
    dem[dem.mask] = np.nan
    create_raster(
        f"{qgis_path}/dem.tiff",
        dem,
        (ORIGIN_X, ORIGIN_Y),
        CRS_EPSG,
        PIXEL_WIDTH,
        PIXEL_HEIGHT,
    )

    # river
    river = extract_post_gauges.getRaster(RIVER_FILE, no_data=-1)
    river_mask = np.zeros(river.shape)
    river_mask[river == -1] = 1
    river = ma.masked_array(river, mask=river_mask)
    river[river.mask] = -1
    create_raster(
        f"{qgis_path}/river.tiff",
        river,
        (ORIGIN_X, ORIGIN_Y),
        CRS_EPSG,
        PIXEL_WIDTH,
        PIXEL_HEIGHT,
        nan_value=-1,
        rdtype=gdal.GDT_Int32,
    )

    landuse = ma.masked_array(
        extract_post_gauges.getRaster(LANDUSE_FILE, -1),
        mask=getCatchmentMask(CATCHMENT_FILE),
    )

    # landuse
    landuse[landuse.mask] = -1
    create_raster(
        f"{qgis_path}/landuse.tiff",
        landuse,
        (ORIGIN_X, ORIGIN_Y),
        CRS_EPSG,
        PIXEL_WIDTH,
        PIXEL_HEIGHT,
        nan_value=-1,
        rdtype=gdal.GDT_Int32,
    )

    # legend for landuse
    level1_colors = {
        # 1: "black",
        1: (0, 0, 0),
        # 2: "yellow",
        2: (255, 255, 0),
        # 3: "green",
        3: (0, 255 / 2, 0),
        # 4: "mediumaquamarine",
        # 5: "blue",
    }

    unique_ids = np.sort(np.unique(landuse[~landuse.mask]))[1:]
    cathegories = {}
    for i in unique_ids:
        cathegory = i // 100
        if cathegory not in cathegories:
            cathegories[cathegory] = []
        cathegories[cathegory].append(i)

    palette = []
    for cathegory, values in cathegories.items():
        palette += list(
            generateColorMap(
                len(values),
                level1_colors[cathegory][0],
                level1_colors[cathegory][1],
                level1_colors[cathegory][2],
                exclude_white=True,
            )
        )

    with open(f"{qgis_path}/landuse_legend.txt", "w") as file:
        for code, color in zip(unique_ids, palette):
            file.write(
                f"{code} {color[0]} {color[1]} {color[2]} 255 {getNameForLanduse(code)}\n"
            )

    # ldd map
    ldd = ma.masked_array(
        extract_post_gauges.getRaster(LDD_FILE, -1),
        mask=getCatchmentMask(CATCHMENT_FILE),
    )
    ldd[ldd.mask] = 255
    create_raster(
        f"{qgis_path}/ldd.tiff",
        ldd,
        (ORIGIN_X, ORIGIN_Y),
        CRS_EPSG,
        PIXEL_WIDTH,
        PIXEL_HEIGHT,
        nan_value=255,
        rdtype=gdal.GDT_Byte,
    )

    # ldd legend
    ldd_legends = {
        1: "Drainage direction to the southwest",
        2: "Drainage direction to the south",
        3: "Drainage direction to the southeast",
        4: "Drainage direction to the west",
        5: "No drainage direction (e.g. a pit)",
        6: "Drainage direction to the east",
        7: "Drainage direction to the northwest",
        8: "Drainage direction to the north",
        9: "Drainage direction to the northeast",
    }
    palette = list(
        generateColorMap(
            9,
            0,
            0,
            0,
            True,
        )
    )
    with open(f"{qgis_path}/ldd_legend.txt", "w") as file:
        for direction, color in zip(range(1, 10), palette):
            file.write(
                f"{direction} {color[0]} {color[1]} {color[2]} 255 {ldd_legends[direction]}\n"
            )

    # streamorder
    streamorder = ma.masked_array(
        extract_post_gauges.getRaster(STREAMORDER_FILE),
        mask=getCatchmentMask(CATCHMENT_FILE),
    )
    streamorder[streamorder.mask] = np.nan
    create_raster(
        f"{qgis_path}/streamorder.tiff",
        streamorder,
        (ORIGIN_X, ORIGIN_Y),
        CRS_EPSG,
        PIXEL_WIDTH,
        PIXEL_HEIGHT,
    )

    # stream order legend
    unique_ids = np.sort(np.unique(streamorder[~streamorder.mask]))[:-1]
    palette = list(
        generateColorMap(
            len(unique_ids),
            0,
            0,
            0,
            True,
        )
    )
    palette.reverse()
    with open(f"{qgis_path}/streamorder_legend.txt", "w") as file:
        for streamorder_id, color in zip(unique_ids, palette):
            file.write(
                f"{streamorder_id:g} {color[0]} {color[1]} {color[2]} 255 Streamorder {streamorder_id:g}\n"
            )


def main() -> None:
    plt.rcParams.update(
        {
            "pgf.preamble": r"\usepackage{siunitx}",
            "font.family": "serif",
            "pgf.rcfonts": False,
        }
    )

    # writeDem(DEM_FILE, CATCHMENT_FILE, DEM_OUTFILE)
    writeRiverAndGauges(DEM_FILE, CATCHMENT_FILE, RIVER_FILE, GAUGE_FILE, RIVER_OUTFILE)
    # writeStreamOrder(CATCHMENT_FILE, STREAMORDER_FILE, STREAMORDER_OUTFILE)
    # writeLdd(CATCHMENT_FILE, LDD_FILE, LDD_OUTFILE)
    # writeLandUse(CATCHMENT_FILE, LANDUSE_FILE, LANDUSE_OUTFILE, LANDUSE_LEVEL_1_OUTFILE)
    prepareQgis()


if __name__ == "__main__":
    main()
