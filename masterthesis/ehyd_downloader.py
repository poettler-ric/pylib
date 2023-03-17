#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:07:09 2023

@author: gegese
"""

from argparse import ArgumentParser
from logging import info
from os import walk
from os.path import join as pjoin
from re import compile
from requests import get as rget
import logging

__ID_DUMP_FILE = "dumped_id.txt"


def try_download(url, out_file):
    r = rget(url, allow_redirects=True)
    if len(r.content) > 0:
        with open(out_file, "wb") as f:
            info(f"writing {url} to file")
            f.write(r.content)


def determine_start_id(basepath, initial_start):
    max_id = initial_start
    csv_pattern = compile("(\\d+)_runoff.csv")
    for _, _, filenames in walk(basepath):
        for filename in filenames:
            match = csv_pattern.match(filename)
            if match:
                max_id = max(max_id, int(match.group(1)))
            if filename == __ID_DUMP_FILE:
                with open(pjoin(basepath, filename)) as f:
                    max_id = max(max_id, int(f.read()))
    return max_id


def brute_downloader(basepath):
    start_id = determine_start_id(basepath, 200_000)
    info(f"starting at id {start_id}")
    end_id = 400_000
    for i in range(start_id, end_id):
        if not i % 500:
            info(f"trying id {i}")
            with open(pjoin(basepath, __ID_DUMP_FILE), "w") as f:
                f.write(str(i))

        try_download(
            f"https://ehyd.gv.at/eHYD/MessstellenExtraData/owf?id={i}&file=4",
            pjoin(basepath, f"{i}_runoff.csv"),
        )


def main():
    logging.getLogger().setLevel(logging.INFO)

    parser = ArgumentParser()
    parser.add_argument("basepath", help="basepath to safe the files to")
    args = parser.parse_args()

    brute_downloader(args.basepath)


if __name__ == "__main__":
    main()
