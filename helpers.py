#!/usr/bin/env python3

__author__ = "Richard Pöttler"
__copyright__ = "Copyright (c) 2020 Richard Pöttler"
__license__ = "MIT"
__email__ = "richard.poettler@gmail.com"


def h1(h):
    print(f"=== {h} ===")


def h2(h):
    print(f"--- {h} ---")


def print_df(h, df, filename=None, quiet=False):
    h1(h)
    if not quiet:
        print(df)
    if filename:
        df.to_csv(filename, index=False)
