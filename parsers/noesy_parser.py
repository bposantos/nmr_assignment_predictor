# io/noesy_parser.py

from .spectrum_objects import Peak


def parse_noesy(filepath):

    peaks = []

    with open(filepath) as f:

        next(f)

        for line in f:

            cols = line.split()

            peaks.append(
                Peak(
                    assignment=cols[0],
                    w1=float(cols[1]),
                    w2=float(cols[2]),
                    intensity=float(cols[3]),
                    spectrum_type="NOESY"
                )
            )

    return peaks