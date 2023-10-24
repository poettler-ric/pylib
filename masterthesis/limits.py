#!/usr/bin/env python3

import statistics
from dataclasses import dataclass, field
from os.path import join as pjoin

import prettytable


@dataclass
class LatexData(object):
    name: str
    unit: str


__LATEX_NAME = {
    "BetaSeepage": LatexData("\\beta", "\\si{-}"),
    "Cflux": LatexData("Cflux", "\\si{\\milli\\meter\\per\\day}"),
    "Cfmax": LatexData("Cfmax", "\\si{\\milli\\meter\\per\\celsius\\day}"),
    "CFR": LatexData("CFR", "\\si{-}"),
    "FC": LatexData("FC", "\\si{\\milli\\meter}"),
    "ICF1": LatexData("ICF", "\\si{\\milli\\meter}"),
    "K0": LatexData("K_{0}", "\\si{\\per\\day}"),
    "K4": LatexData("K_{4}", "\\si{\\per\\day}"),
    "KQuickFlow": LatexData("K_{QuickFlow}", "\\si{\\per\\day}"),
    "LP": LatexData("LP", "\\si{-}"),
    "N": LatexData("N", "\\si{\\second\\per\\sqrt[3]{\\m}}"),
    "PERC": LatexData("PERC", "\\si{\\milli\\meter\\per\\day}"),
    "SUZ": LatexData("S_{UZ}", "\\si{\\milli\\meter}"),
    "TT": LatexData("TT", "\\si{\\celsius}"),
    "TTI": LatexData("TTI", "\\si{\\celsius}"),
    "WHC": LatexData("WHC", "\\si{\\milli\\meter\\per\\milli\\meter}"),
}


@dataclass
class Limit(object):
    name: str
    start: float
    end: float
    calibrated: dict[str, float] = field(default_factory=dict)


@staticmethod
def parse_limits(
    lookup_definiton: dict[str, str], inmap_folder: str
) -> dict[str, Limit]:
    limits = {}
    for key, value in lookup_definiton.items():
        factors = value.split(">")
        from_factor, to_factor = float(factors[0]), float(factors[1])

        with open(pjoin(inmap_folder, f"{key}.tbl"), encoding="utf-8") as tbl_file:
            lines = tbl_file.readlines()
            if len(lines) > 1:
                values = []
                for line in lines:
                    values.append(float(line.split()[-1]))
                limits[key] = Limit(
                    key, min(values) * from_factor, max(values) * to_factor
                )
            else:
                line = lines.pop()
                center = float(line.split()[-1])
                limits[key] = Limit(key, center * from_factor, center * to_factor)
    return limits


def parse_calibrated(
    limits: dict[str, Limit], name: str, intbl: str
) -> dict[str, Limit]:
    for key, value in limits.items():
        result = 0.0
        with open(pjoin(intbl, f"{key}.tbl"), encoding="utf-8") as tbl_file:
            lines = tbl_file.readlines()
            if len(lines) > 1:
                values = []
                for line in lines:
                    values.append(float(line.split()[-1]))
                result = statistics.mean(values)
            else:
                line = lines.pop()
                result = float(line.split()[-1])
            pass
        value.calibrated[name] = result
    return limits


@staticmethod
def print_limits_latex(limits: dict[str, Limit]) -> None:
    for _, limit in limits.items():
        print(
            f"${__LATEX_NAME[limit.name].name}$ & {__LATEX_NAME[limit.name].unit} & {limit.start:3f} & {limit.end:3g} \\\\"
        )


@staticmethod
def print_limits(limits: dict[str, Limit]) -> None:
    table = prettytable.PrettyTable()
    table.float_format = ".3"
    table.align["From"] = "r"
    table.align["To"] = "r"
    table.field_names = ["Name", "From", "To"]
    for _, limit in limits.items():
        table.add_row([limit.name, limit.start, limit.end])
    print(table)


def main():
    opt_lookup = {
        # doesn't matter, since SetKquickFlow = 1
        # "AlphaNL": "0.1>6.0",
        "BetaSeepage": "0.01>10.0",
        "Cflux": "0.0>3.0",
        "Cfmax": "0.1>3.0",
        # TODO: reactivate
        # "CFR": "0.0>0.1",
        "FC": "0.1>5.0",
        # doesn't matter, since SetKquickFlow = 1
        # "HQ": "0.01>5.0",
        "ICF1": "0.1>20",
        "K0": "0.2>13.0",
        "K4": "0.05>10.0",
        # doesn't matter, since SetKquickFlow = 1
        # "KHQ": "0.5>8.0",
        "KQuickFlow": "0.5>8.0",
        "LP": "0.1>1.8",
        # nowhere found in code
        # "MaxLeakage": "0.0>3.0",
        "N": "0.5>2.0",
        "PERC": "0.1>50.0",
        "SUZ": "0.1>1.0",
        "TT": "-2.0>2.0",
        "TTI": "0.0>7.0",
        "WHC": "0.1>3.0",
    }

    calibration_intbl = "/data/home/richi/master_thesis/model_MAR/intbl"
    intbl_inca = "/data/home/richi/master_thesis/model_MAR/intbl_calib_INCA/"
    intbl_era5_nse = "/data/home/richi/master_thesis/model_MAR/intbl_calib_ERA5_NSE"
    intbl_era5_kge = "/data/home/richi/master_thesis/model_MAR/intbl_calib_ERA5_KGE"

    limits = parse_limits(opt_lookup, calibration_intbl)
    parse_calibrated(limits, "inca_nse", intbl_inca)
    parse_calibrated(limits, "era5_nse", intbl_era5_nse)
    parse_calibrated(limits, "era5_kge", intbl_era5_kge)

    print_limits(limits)
    print_limits_latex(limits)


if __name__ == "__main__":
    main()
