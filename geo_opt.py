#!/usr/bin/env python

# %run qmmm_ONIOM_code.py -fstrc nC18-rs/structure_num_1_of_nc18_TS_ONIOM_rs_dd.pdb -reff1 nC18-rs/nC18_rs_PBE_cc_UFF.com -reff2 nC18-rs/nC18_rs_PBE_aug_cc_UFF.com -tol1 0.08 -tol2 0.05 -gl1 nC18-rs/gaussian1 -gl2 nC18-rs/gaussian2 -xyz_fl1 nC18-rs/traj1_opt.xyz -xyz_fl2 nC18-rs/traj2_opt.xyz -d nC18-rs/ -rxn_type rs -fl1 nC18-rs/stdout1.log -fl2 nC18-rs/stdout2.log -ftr1 nC18-rs/stdout1.traj -ftr2 nC18-rs/stdout2.traj

#!/usr/bin/env python

# Gaussian gives force_consistent parameter as False

# things still needing work
# 1) how to input memory information
# 2) freeze atoms (make forces zero on frozen atoms) (done)


"""Use ASE (git: fea26c) to perform geometry optimization using damped dynamics."""
# Need two step optimization
# STEP 1 - pVDZ, fmax < 0.08
# STEP 2 - AUG-cc-pVDZ, fmax < 0.05

from __future__ import print_function, unicode_literals
from builtins import str, range
import argparse

# import mpi4py # including this gives mpi error
from ase.io import read, Trajectory, write, trajectory
import ase.parallel as mpi
import numpy as np
from ase import Atoms
from ase.calculators.gaussian import Gaussian
import numpy as np
from ase.optimize import FIRE as DampedDynamics
import Gaussian_ONIOM_plus_ASE_for_QMMM.modules.gaussian_ONIOM_calculator as goc
from ase.io import read


def get_gaussian_params():  # this function is not needed for calculation but needed for setup (to create the calculator)
    params = dict(
        chk="gaussian.chk",
        # xc='PBEPBE',
        basis="Gen/Auto",
        # basisfile='gen.basis',
        # SCRF='SMD,Solvent=Water',
        # Geom='(Check, NewDefinition)',
        # Guess='Read',
        Integral="Grid=UltraFine",
        SCF="Tight,IntRep,XQC",
        Symmetry="None",
        pop="Hirshfeld",
        EmpiricalDispersion="GD3BJ",
    )

    # if 'nprocshared' not in params:
    #    params['nprocshared'] = sys.core_count

    # if 'mem' not in params and sys.mem_bytes is not None:
    #    params['mem'] = '{}GB'.format(int(sys.mem_bytes*0.6/(1024.**3)/(sys.core_count/params['nprocshared'])))

    print(params)
    return params


def get_calculator(calc_params, filename, calc_label, mem, nproc):
    return goc.Gaussian_ONIOM(
        memory=mem, nprocshar=nproc, fname=filename, label=calc_label, **calc_params
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Geometry optimization using Gaussian ONIOM"
    )
    parser.add_argument(
        "-fstrc", "--fname_strc", action="store", required=True, default="abc.pdb"
    )  # input structure to perform optimization on
    parser.add_argument(
        "-reff1", "--ref_file1", action="store", required=True, default="Gaussian.com"
    )  # ref1 gaussian_file made by gaussview
    parser.add_argument(
        "-reff2", "--ref_file2", action="store", required=True, default="Gaussian.com"
    )  # ref2 gaussian_file made by gaussview
    #    parser.add_argument('-f', '--fname', action='store', required=True, default="gaussian") # prefix of gaussian input_file
    parser.add_argument(
        "-d", "--dir", action="store", required=True
    )  # directory to store STEP1 opt_file
    parser.add_argument(
        "-gl1", "--gaus_label1", action="store", type=str, default="gaussian1"
    )  # STEP 1 gaussian output file name
    parser.add_argument(
        "-gl2", "--gaus_label2", action="store", type=str, default="gaussian2"
    )  # STEP 2 gaussian output file name
    parser.add_argument(
        "-tol1", "--tolerance1", action="store", required=True, default=0.08, type=float
    )  # force convergence criteria tolerance for STEP 1
    parser.add_argument(
        "-tol2", "--tolerance2", action="store", required=True, default=0.05, type=float
    )  # force convergence criteria tolerance for STEP 2
    parser.add_argument(
        "-xyz_fl1",
        "--xyz_file_name1",
        action="store",
        required=True,
        default="traj1.xyz",
    )  # file to store trajectory1 in xyz format
    parser.add_argument(
        "-xyz_fl2",
        "--xyz_file_name2",
        action="store",
        required=True,
        default="traj2.xyz",
    )  # file to store trajectory1 in xyz format
    parser.add_argument(
        "-fl1", "--log_filename1", action="store", default="stdout1"
    )  # log file for STEP1
    parser.add_argument(
        "-fl2", "--log_filename2", action="store", default="stdout2"
    )  # log file for STEP2
    parser.add_argument(
        "-ftr1", "--traj_filename1", action="store", default="traj1"
    )  # traj file for STEP1
    parser.add_argument(
        "-ftr2", "--traj_filename2", action="store", default="traj2"
    )  # traj file for STEP2
    # parser.add_argument('-rxn_type', '--reaction_type', action='store', required=True, default='rs') # rs for reactant state and ps for product state
    parser.add_argument("-mem", "--memory", action="store", default=1, type=int)
    parser.add_argument("-nproc", "--nprocshared", action="store", default=1, type=int)
    args = parser.parse_args()
    return args


def run_optimizer(dd_optimizer, tol):
    #    dd = DampedDynamics(dd_optimizer)
    dd_optimizer.run(fmax=tol, steps=1000)
    return dd_optimizer


def set_optimizer(atoms, filename, traj_file_name):
    dd_optimizer = DampedDynamics(
        atoms, restart=None, logfile=filename, trajectory=traj_file_name
    )
    return dd_optimizer


if __name__ == "__main__":
    args = parse_args()
    if mpi.rank == 0:
        print(args)

    # Calculator parameters
    # calc_params = get_gaussian_params(
    #    charge=args.charge,
    #    multiplicity=args.multiplicity,
    #    nprocshared=args.proc_per_image,
    # )

    # calc_label1 = args.gaus_label1

    # if mpi.size > 1:
    #    calc_label = str(mpi.rank) + '-' + calc_label
    # if len(args.tmp) > 0:
    #    calc_label = args.tmp.rstrip('/') + '/' + calc_label

    calc_params = get_gaussian_params()

    atoms1 = read(args.fname_strc)
    atoms1.calc = get_calculator(
        calc_params, args.ref_file1, args.gaus_label1, args.memory, args.nprocshared
    )
    # STEP 1

    dd_optimizer1 = set_optimizer(atoms1, args.log_filename1, args.traj_filename1)
    dd_optimizer1 = run_optimizer(dd_optimizer1, tol=args.tolerance1)
    ########

    p = read(args.traj_filename1)
    STEP1_opt_fl_name = args.dir + "STEP1_opt.pdb"  # set it as user input
    p.write(STEP1_opt_fl_name)
    # atom_types_ret = atoms1.arrays["atomtypes"]

    # STEP2
    atoms2 = read(STEP1_opt_fl_name)
    atoms2.calc = get_calculator(
        calc_params, args.ref_file2, args.gaus_label2, args.memory, args.nprocshared
    )
    dd_optimizer2 = set_optimizer(atoms2, args.log_filename2, args.traj_filename2)
    dd_optimizer2 = run_optimizer(dd_optimizer2, tol=args.tolerance2)
    ########

    p = read(args.traj_filename2)
    STEP2_opt_fl_name = args.dir + "STEP2_opt.pdb"  # set it as user input
    p.write(STEP2_opt_fl_name)

    # to convert .traj file to .xyz file

    k = trajectory.TrajectoryReader(args.traj_filename1)
    write(args.xyz_file_name1, k)

    k = trajectory.TrajectoryReader(args.traj_filename2)
    write(args.xyz_file_name2, k)
