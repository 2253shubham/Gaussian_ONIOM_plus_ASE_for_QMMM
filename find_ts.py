#!/usr/bin/env python

"""Use ASE (git: fea26c) to perform a transition-state search."""


from __future__ import print_function, unicode_literals
from builtins import str, range
import argparse
import mpi4py
from ase.calculators.gaussian import Gaussian
from ase.io import read, Trajectory, write, trajectory
from ase.neb import NEB
from ase.optimize import FIRE as DampedDynamics
from ase.optimize import BFGS as QuasiNewton
import ase.parallel as mpi
import Gaussian_ONIOM_plus_ASE_for_QMMM.modules.gaussian_ONIOM_calculator as goc
import numpy as np
# from asemc.util import sys


def get_atoms_from_file(filename, calc_params, calc_label):
    atoms = read(filename)
    atoms.set_calculator(get_calculator(calc_params, calc_label))
    atoms.cell = 3e1 * np.eye(3, dtype=np.float_)
    # atoms.pbc = (True, True, True)
    return atoms


def get_calculator(calc_params, filename, calc_label, mem, nproc):
    return goc.Gaussian_ONIOM(
        memory=mem, nprocshar=nproc, fname=filename, label=calc_label, **calc_params
    )


def get_gaussian_params(**kwargs):
    params = dict(
        chk="gaussian.chk",
        method="UM06L",
        basis="Gen/Auto",
        basisfile="gen.basis",
        SCRF="SMD,Solvent=Water",
        # Geom='(Check, NewDefinition)',
        # Guess='Read',
        Integral="(Grid=UltraFine)",
        SCF="(Tight,IntRep,XQC)",
        Symmetry=None,
        **kwargs
    )

    if "nprocshared" not in params:
        params["nprocshared"] = sys.core_count

    if "mem" not in params and sys.mem_bytes is not None:
        params["mem"] = "{}GB".format(
            int(
                sys.mem_bytes
                * 0.6
                / (1024.0**3)
                / (sys.core_count / params["nprocshared"])
            )
        )

    return params


def neb_data(path, distance, energies, fname):
    """Write out reaction coordinates and energies."""

    data = []
    for i, apath in enumerate(path):
        data.append([np.sqrt(distance(apath, path[0])[0]), energies[i] - energies[0]])
    np.savetxt(fname, np.asarray(data), delimiter=" ")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Automatic transition state" " search (ReaxFF/Gaussian)."
    )
    parser.add_argument(
        "-reff1",
        "--ref_gau_file1",
        action="store",
        required=True,
        default="Gaussian.com",
    )  # ref1 gaussian_file made by gaussview
    parser.add_argument(
        "-reff2",
        "--ref_gau_file2",
        action="store",
        required=True,
        default="Gaussian.com",
    )  # ref2 gaussian_file made by gaussview
    parser.add_argument(
        "-gl1", "--gaus_label1", action="store", type=str, default="gaussian1"
    )  # STEP 1 gaussian output file name
    parser.add_argument(
        "-gl2", "--gaus_label2", action="store", type=str, default="gaussian2"
    )  # STEP 2 gaussian output file name
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-ni", "--neb_l", action="store", default=-1, type=int)
    parser.add_argument(
        "-ngc", "--neb_guide_count", action="store", default=0, type=int
    )
    parser.add_argument("-ngp", "--neb_guide_prefix", action="store", default="neb")
    parser.add_argument("-ngs", "--neb_guide_suffix", action="store", default=".pdb")
    parser.add_argument("-ps", "--prod_struct", action="store")
    parser.add_argument("-ppi", "--proc-per-image", action="store", default=1, type=int)
    parser.add_argument("-rs", "--react_struct", action="store", required=True)
    parser.add_argument("-mem", "--memory", action="store", default=1, type=int)
    parser.add_argument("-nproc", "--nprocshared", action="store", default=1, type=int)
    parser.add_argument("--tmp", action="store", default="")
    args = parser.parse_args()
    return args


def prepare_neb(
    atoms1,
    atoms2,
    ref_file,
    calc_label,
    mem,
    nproc,
    nimages,
    calc_params,  # calc_label,
    nguide=0,
    prefix="neb",
    suffix=".pdb",
):
    images = [atoms1]
    guide_images = [0]

    for i in range(nimages):
        if i < nguide:
            guide_images.append(i + 1)
            int_im = read(prefix + str(i + 1) + suffix)
            # int_im.mem = calc_params['mem']
            # int_im.nproc = calc_params['nprocshared']
            int_im.calc = get_calculator(
                calc_params, ref_file, calc_label.format(id=i + 1), mem, nproc
            )
            images.append(int_im)
        else:
            at = atoms1.copy()
            # at.mem = calc_params['mem']
            # at.nproc = calc_params['nprocshared']
            at.calc = get_calculator(
                calc_params, ref_file, calc_label.format(id=i + 1), mem, nproc
            )
            images.append(at)

    images.append(atoms2)
    guide_images.append(-1)

    neb = NEB(
        images,
        k=0.1,
        climb=False,
        method="improvedtangent",
        remove_rotation_and_translation=False,
        parallel=True,
    )
    neb.interpolate(guide_images=guide_images)

    if mpi.rank == 0:
        write_neb(images)

    return neb


def run_neb(neb, tol=0.08):
    tol_stage1 = tol

    dd = DampedDynamics(neb)
    qn = QuasiNewton(neb)

    for i, im in enumerate(neb.images[1:-1]):
        if i == mpi.rank:
            traj = Trajectory("neb{}.traj".format(i + 1), "w", im, master=True)
            dd.attach(traj)
            qn.attach(traj)

    dd.run(fmax=tol_stage1, steps=1)

    if mpi.rank == 0:
        print('*Add "Guess=Read"')

    for im in neb.images[1:-1]:
        im.calc.set(Guess="Read")

    dd.run(fmax=tol_stage1, steps=1000)

    # if mpi.rank == 0:
    #     print('== Switch to QuasiNewton (BFGS) ==')

    # qn.run(fmax=tol, steps=500)

    return neb


def write_neb(images=None):
    """Convert an ASE trajectory to an ARC movie file."""

    if images is None:
        # images = read('neb.traj@-18:')
        images = read("neb.traj@:18")

    from asemc.io.helper import write_movie

    # write_movie('neb.arc', path, xyz.atomtypes,
    #             SimulationSystem.calc.get_boxvec(), shift=[0, 0, 0.5])
    write_movie(
        "neb.res",
        [x.positions for x in images],
        [x.get_chemical_symbols() for x in images],
        [x.cell for x in images],
        filetype="res",
        shift=[0.5, 0.5, 0.5],
    )


def write_movie_alt(filename, nimages, index):
    # open('latest_update.pdb', 'w').close()
    for i in range(nimages):
        tr = trajectory.TrajectoryReader("neb" + str(i + 1) + ".traj")
        write(filename, tr[index], append=True)


if __name__ == "__main__":
    args = parse_args()
    if mpi.rank == 0:
        print(args)

    # Calculator parameters
    calc_params = get_gaussian_params(
        charge=0,
        multiplicity=1,
        nprocshared=args.proc_per_image,
    )

    calc_label = "{id}/gaussian"
    if mpi.size > 1:
        calc_label = str(mpi.rank) + "-" + calc_label
    if len(args.tmp) > 0:
        calc_label = args.tmp.rstrip("/") + "/" + calc_label

    # Create the Atoms object
    min1 = read(args.react_struct)
    # min1.mem = calc_params['mem']
    # min1.nproc = calc_params['nprocshared']
    min1.calc = get_calculator(
        calc_params,
        args.ref_gau_file1,
        calc_label.format(id=0),
        args.memory,
        args.nprocshared,
    )

    if args.prod_struct is not None:
        min2 = read(args.prod_struct)
        # min2.mem = calc_params['mem']
        # min2.nproc = calc_params['nprocshared']
        min2.calc = get_calculator(
            calc_params,
            args.ref_gau_file1,
            calc_label.format(id=args.neb_l + 1),
            args.memory,
            args.nprocshared,
        )

        if args.neb_l > 2:
            neb = prepare_neb(
                min1,
                min2,
                args.ref_gau_file1,
                calc_label,
                args.memory,
                args.nprocshared,
                nimages=args.neb_l,
                calc_params=calc_params,
                #  calc_label=calc_label,
                nguide=args.neb_guide_count,
                prefix=args.neb_guide_prefix,
                suffix=args.neb_guide_suffix,
            )
            if not args.dry_run:
                if mpi.size != args.neb_l:
                    raise RuntimeError("# of MPI ranks != # of NEB images")
                neb = run_neb(neb)

    k = trajectory.TrajectoryReader("neb" + str(3) + ".traj")
    for i in range(len(k)):
        write_movie_alt("neb_step_" + str(i + 1) + "_movie.pdb", args.neb_l, i)
