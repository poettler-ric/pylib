#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 14:01:17 2023

@author: gegese
"""

from argparse import ArgumentParser
import fiona
import pyproj
import os


def get_metadata(text, coordrefs):
    include = True

    meta = {
        "state": "",
        "ID": 0,
        "name": "",
        "area": 0,
        "elevation": 0,
        "x": 0,
        "y": 0,
        "unit": "",
    }

    end_found = False

    i = 0
    while not end_found:
        # print(text[i])

        # check state
        if "HZB-Nummer:" in text[i]:
            line = int(text[i].split(";")[1])
            meta["ID"] = line

        if "Messstelle:" in text[i]:
            line = text[i].split(";")[1]
            meta["name"] = line

        # check state
        if "Dienststelle:" in text[i]:
            line = text[i].split(";")[1]
            coordrefs[line]
            meta["state"] = line

        # check for area
        if "orogr.Einzugsgebiet" in text[i]:
            line = float(text[i].split(";")[1].replace(",", "."))
            meta["area"] = line

        # check for elevation
        j = 0
        if "Pegelnullpunkt" in text[i]:
            end_gauge = False
            while not end_gauge:
                if "Bundesmeldenetz" in text[i + j]:
                    end_gauge = True
                    # text
                    line = float(
                        text[i + j - 1].split(";")[1].replace(",", ".").lstrip("0")
                    )
                    meta["elevation"] = line
                j = j + 1

        # check for x and y
        j = 0
        if "Koordinaten" in text[i]:
            end_gauge = False
            while not end_gauge:
                if "Exportzeitreihe" in text[i + j]:
                    end_gauge = True
                    # text
                    line = text[i + j - 1].split(";")[1].split("-")
                    x = float(line[0].lstrip(" ").rstrip(" ").replace(",", "."))
                    y = float(line[1].lstrip(" ").rstrip(" ").replace(",", "."))
                    meta["x"] = x
                    meta["y"] = y
                j = j + 1

        # check if end is reached
        if "Werte:" in text[i]:
            end_found = True

        # check for unit
        if "Einheit:" in text[i]:
            line = text[i].split(";")[1]
            if "m" and "/" and "s" in line:
                meta["unit"] = "discharge"
            else:
                meta["unit"] = "false"
                include = False

        i += 1

    return meta, include


def main():
    parser = ArgumentParser()
    parser.add_argument("directory", help="location of the metadata files")
    parser.add_argument("shapefile", help="location of the resulting shape file")
    args = parser.parse_args()

    # reference systems
    coordrefs = {
        "HD-Vorarlberg": "+proj=tmerc +lat_0=0 +lon_0=10.3333333333333 +k=1 +x_0=150000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        "HD-Tirol": "+proj=tmerc +lat_0=0 +lon_0=10.3333333333333 +k=1 +x_0=150000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        "HD-Salzburg": "+proj=tmerc +lat_0=0 +lon_0=13.3333333333333 +k=1 +x_0=450000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        r"HD-Oberösterreich": "+proj=tmerc +lat_0=0 +lon_0=13.3333333333333 +k=1 +x_0=450000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        r"HD-Kärnten": "+proj=tmerc +lat_0=0 +lon_0=13.3333333333333 +k=1 +x_0=450000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        r"HD-Niederösterreich": "+proj=tmerc +lat_0=0 +lon_0=16.3333333333333 +k=1 +x_0=750000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        "HD-Burgenland": "+proj=tmerc +lat_0=0 +lon_0=16.3333333333333 +k=1 +x_0=750000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        "HD-Steiermark": "+proj=tmerc +lat_0=0 +lon_0=16.3333333333333 +k=1 +x_0=750000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        "via donau (WSD)": "+proj=tmerc +lat_0=0 +lon_0=16.3333333333333 +k=1 +x_0=750000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        "HD-Wien (MA 45)": "+proj=tmerc +lat_0=0 +lon_0=16.3333333333333 +k=1 +x_0=750000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
        "Magistratsabteilung 31": "+proj=tmerc +lat_0=0 +lon_0=16.3333333333333 +k=1 +x_0=750000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs",
    }

    # target crs
    target_crs = "+proj=lcc +lat_0=47.5 +lon_0=13.3333333333333 +lat_1=49 +lat_2=46 +x_0=400000 +y_0=400000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs"
    # shape settings
    schema = {
        "geometry": "Point",
        "properties": [
            ("state", "str"),
            ("ID", "int"),
            ("name", "str"),
            ("area", "float"),
            ("elevation", "float"),
        ],
    }
    pointShp = fiona.open(
        args.shapefile,
        mode="w",
        driver="ESRI Shapefile",
        schema=schema,
        crs="EPSG:31287",
    )

    # Main loop
    for file in os.listdir(args.directory):
        filename = os.fsdecode(file)

        print(filename)

        if "runoff" in filename:
            # construct path
            path = args.directory + "/" + filename
            # open file
            f = open(path, "r")
            text = f.read().splitlines()
            # get metadata
            meta, include = get_metadata(text, coordrefs)
            # perform projection
            source_crs = coordrefs[meta["state"]]

            # special case for east tirol
            if meta["state"] == "HD-Tirol" and meta["x"] > 275000:
                source_crs = "+proj=tmerc +lat_0=0 +lon_0=13.3333333333333 +k=1 +x_0=450000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs"

            if meta["ID"] == 201830:
                source_crs = "+proj=tmerc +lat_0=0 +lon_0=10.3333333333333 +k=1 +x_0=150000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs"

            if meta["ID"] == 212217:
                source_crs = "+proj=tmerc +lat_0=0 +lon_0=10.3333333333333 +k=1 +x_0=150000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs"

            if meta["ID"] == 207035:
                source_crs = "+proj=tmerc +lat_0=0 +lon_0=13.3333333333333 +k=1 +x_0=450000 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232 +units=m +no_defs"

            transfo = pyproj.Transformer.from_crs(source_crs, target_crs)
            x_new, y_new = transfo.transform(meta["x"], meta["y"])
            # check if the point should be included
            if include:
                rowDict = {
                    "geometry": {"type": "Point", "coordinates": (x_new, y_new)},
                    "properties": {
                        "state": meta["state"],
                        "ID": meta["ID"],
                        "name": meta["name"],
                        "area": meta["area"],
                        "elevation": meta["elevation"],
                    },
                }
                pointShp.write(rowDict)

    pointShp.close()


if __name__ == "__main__":
    main()
