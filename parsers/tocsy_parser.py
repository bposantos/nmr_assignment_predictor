# io/tocsy_parser.py

from .spectrum_objects import Peak


def parse_tocsy(filepath):

    peaks = []

    with open(filepath) as f:

        next(f)

        for line in f:

            cols = line.split()

            assignment = cols[0]
            w1 = float(cols[1])
            w2 = float(cols[2])
            intensity = float(cols[3])

            peaks.append(
                Peak(
                    assignment=assignment,
                    w1=w1,
                    w2=w2,
                    intensity=intensity,
                    spectrum_type="TOCSY"
                )
            )

    return peaks
