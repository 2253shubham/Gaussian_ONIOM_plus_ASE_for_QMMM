#!/usr/bin/env python

"""Use ASE (git: fea26c) to perform a transition-state search."""


from __future__ import print_function, unicode_literals
from builtins import str, range
import argparse
from ase.io import trajectory, write, read
import ase.parallel as mpi
import numpy as np


def get_atoms_from_file(filename):
    atoms = trajectory.TrajectoryReader(filename)
    return atoms


def parse_args():
    parser = argparse.ArgumentParser(description="Creating intermediate pdb files.")
    parser.add_argument(
        "-d", "--dir", action="store", required=True
    )  # directory to store neb.pdb files
    parser.add_argument(
        "-n", "--num_itd", action="store", default=8, type=int, required=True
    )
    parser.add_argument(
        "-traj", "--traj_file", action="store", required=True, default="traj.traj"
    )  # trajectory input file
    parser.add_argument(
        "-rxn_type", "--reaction_type", action="store", required=True, default="rs"
    )  # rs for reactant state and ps for product state
    args = parser.parse_args()
    return args


def write_pdb(filename, atoms):
    """Write the atoms object to a pdb file"""

    write(filename, atoms)


if __name__ == "__main__":
    args = parse_args()
    if mpi.rank == 0:
        print(args)

    atoms = get_atoms_from_file(args.traj_file)
    l = len(atoms)
    if args.reaction_type == "rs":
        c = args.num_itd
        for i in range(int(l / args.num_itd), l, int(l / args.num_itd)):
            write_pdb(args.dir + "neb" + str(c) + ".pdb", atoms[i])
            c = c - 1
    else:
        c = args.num_itd + 2
        for i in range(int(l / args.num_itd), l, int(l / args.num_itd)):
            write_pdb(args.dir + "neb" + str(c) + ".pdb", atoms[i])
            c = c + 1
