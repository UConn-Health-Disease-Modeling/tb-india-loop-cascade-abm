from collections import namedtuple
import numpy as np
import os
from argparse import ArgumentParser
from emod_api.spatialreports.spatial import SpatialReport

SCRIPT_PATH = os.path.realpath(__file__)
WORKING_DIRECTORY = os.path.dirname(SCRIPT_PATH)

def main(filename: str):

    report = SpatialReport(filename)

    import matplotlib.pyplot as plt

    y_axis_guess = os.path.basename(filename).strip( "SpatialReport_" ).strip( ".bin" )
    # Show node-wise time series data
    x = [ x for x in range(len(report[1].data)) ]
    # Basic stacked area chart.
    y_values = []
    for node_id in report.node_ids:
    #    plt.plot(report[node_id].data, label=f"Node {node_id}")
        y_values.append( report[node_id].data )
    plt.stackplot(x, y_values )
    plt.title(f"{y_axis_guess} for Nodes")
    plt.legend()
    plt.xlabel("Time Step")
    plt.ylabel(f"{y_axis_guess}")
    plt.show()

    return


def dump_source():

    with open(SCRIPT_PATH, "r") as file:
        print(file.read())

    return


if __name__ == "__main__":
    parser = ArgumentParser() 
    parser.add_argument(
        "-f", "--filename", default=None, help="spatial report filename"
    )
    parser.add_argument(
        "-g",
        "--get",
        default=False,
        action="store_true",
        help="Write source code for this example to stdout.",
    )

    args = parser.parse_args()

    if args.get:
        dump_source()
    elif args.filename:
        main(args.filename)
    else:
        parser.print_help()



