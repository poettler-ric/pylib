#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:07:09 2023

@author: gegese
"""

from argparse import ArgumentParser
from os.path import join as pjoin
from requests import get as rget


def try_download(url, out_file):
    r = rget(url, allow_redirects=True)
    if len(r.content) > 0:
        with open(out_file, "wb") as f:
            f.write(r.content)


def brute_downloader(basepath):
    for i in range(200_000, 400_000):
        try_download(
            f"https://ehyd.gv.at/eHYD/MessstellenExtraData/owf?id={i}&file=1",
            pjoin(basepath, f"{i}_meta.csv"),
        )
        try_download(
            f"https://ehyd.gv.at/eHYD/MessstellenExtraData/owf?id={i}&file=4",
            pjoin(basepath, f"{i}_runoff.csv"),
        )


def main():
    parser = ArgumentParser()
    parser.add_argument("basepath", help="basepath to safe the files to")
    args = parser.parse_args()

    brute_downloader(args.basepath)


if __name__ == "__main__":
    main()
