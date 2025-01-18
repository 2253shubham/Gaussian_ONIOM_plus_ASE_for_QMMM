import numpy as np
from ase.units import Bohr, Hartree
import subprocess
from ase.calculators.gaussian import Gaussian
import re
from ase.calculators.calculator import CalculationFailed
from ase import Atoms
import os
from ase.calculators.singlepoint import SinglePointCalculator
from ase.calculators.calculator import FileIOCalculator, EnvironmentError
from shutil import which
from ase.io.formats import filetype, get_ioformat
from typing import (
    IO, List, Any, Iterable, Tuple, Union, Sequence, Dict, Optional)
import os
import sys 

counter=np.zeros(100) # can be changed based on number of images, to track the number of function calls
np.set_printoptions(threshold=sys.maxsize)

_re_atom = re.compile(
    r'^\s*\S+\s+(\S+)\s+(?:\S+\s+)?(\S+)\s+(\S+)\s+(\S+)\s*$'
)
_re_forceblock = re.compile(r'^\s*Center\s+Atomic\s+Integrated\s+Forces\s+\S+\s*$')
_re_l716 = re.compile(r'^\s*\(Enter .+l716.exe\)$')


def modify_forces(calc_output, frozen):
    raw_forces = calc_output.calc.results["forces"]
    for i in range(len(raw_forces)):
        if frozen[i] != 0:
            raw_forces[i] = [0]*len(raw_forces[i])
    calc_output.calc.results["forces"] = raw_forces


def read_gaussian_out(filename):
    fd = open(filename, "rt")
    atoms = None
    energy = None
    dipole = None
    forces = None
    for line in fd:
        line = line.strip()
        if line.startswith(r'1\1\GINC'):
            # We've reached the "archive" block at the bottom, stop parsing
            break

        if (line == 'Input orientation:'
                or line == 'Standard orientation:'
                or line == 'Z-Matrix orientation:'):

            numbers = []
            positions = []
            pbc = np.zeros(3, dtype=bool)
            cell = np.zeros((3, 3))
            npbc = 0
            # skip 4 irrelevant lines
            for _ in range(4):
                fd.readline()
            while True:
                match = _re_atom.match(fd.readline())
                if match is None:
                    break
                number = int(match.group(1))
                pos = list(map(float, match.group(2, 3, 4)))
                if number == -2:
                    pbc[npbc] = True
                    cell[npbc] = pos
                    npbc += 1
                else:
                    numbers.append(max(number, 0))
                    positions.append(pos)
            atoms = Atoms(numbers, positions, pbc=pbc, cell=cell)
        elif (line.startswith('ONIOM: extrapolated energy =')):
            # Some semi-empirical methods (Huckel, MINDO3, etc.),
            # or SCF methods (HF, DFT, etc.)
            energy = float(line.split('=')[1].split()[0].replace('D', 'e'))
            energy *= Hartree
        elif (line.startswith('E2 =') or line.startswith('E3 =')
                or line.startswith('E4(') or line.startswith('DEMP5 =')
                or line.startswith('E2(')):
            # MP{2,3,4,5} energy
            # also some double hybrid calculations, like B2PLYP
            energy = float(line.split('=')[-1].strip().replace('D', 'e'))
            energy *= Hartree
        elif line.startswith('Wavefunction amplitudes converged. E(Corr)'):
            # "correlated method" energy, e.g. CCSD
            energy = float(line.split('=')[-1].strip().replace('D', 'e'))
            energy *= Hartree
        elif _re_l716.match(line):
            # Sometimes Gaussian will print "Rotating derivatives to
            # standard orientation" after the matched line (which looks like
            # "(Enter /opt/gaussian/g16/l716.exe)", though the exact path
            # depends on where Gaussian is installed). We *skip* the dipole
            # in this case, because it might be rotated relative to the input
            # orientation (and also it is numerically different even if the
            # standard orientation is the same as the input orientation).
            line = fd.readline().strip()
            if not line.startswith('Dipole'):
                continue
            dip = line.split('=')[1].replace('D', 'e')
            tokens = dip.split()
            dipole = []
            # dipole elements can run together, depending on what method was
            # used to calculate them. First see if there is a space between
            # values.
            if len(tokens) == 3:
                dipole = list(map(float, tokens))
            elif len(dip) % 3 == 0:
                # next, check if the number of tokens is divisible by 3
                nchars = len(dip) // 3
                for i in range(3):
                    dipole.append(float(dip[nchars * i:nchars * (i + 1)]))
            else:
                # otherwise, just give up on trying to parse it.
                dipole = None
                continue
            # this dipole moment is printed in atomic units, e-Bohr
            # ASE uses e-Angstrom for dipole moments.
            dipole = np.array(dipole) * Bohr
        elif _re_forceblock.match(line):
            # skip 2 irrelevant lines
            fd.readline()
            fd.readline()
            forces = []
            while True:
                match = _re_atom.match(fd.readline())
                if match is None:
                    break
                forces.append(list(map(float, match.group(2, 3, 4))))
            forces = np.array(forces) * Hartree / Bohr
    if atoms is not None:
        atoms.calc = SinglePointCalculator(
            atoms, energy=energy, dipole=dipole, forces=forces,
        )
    return atoms


class Gaussian_ONIOM(Gaussian):
    implemented_properties = ['energy', 'forces', 'dipole']
    command = 'GAUSSIAN < PREFIX.com > PREFIX.log'
    discard_results_on_any_change = True


    def __init__(self, *args, memory=100, nprocshar=16, fname = 'Gaussian.com', label='Gaussian', **kwargs):
        Gaussian.__init__(self, *args, label=label, **kwargs)
        self.filename_ref = fname
        self.memory = memory
        self.nprocshared = nprocshar


    def write_input(self, atoms, properties=None, system_changes=None):
        label = self.label
        c = 0
        mod_lines = []
        frozen = []
        open_file = open(self.filename_ref, "rt")
        lines = open_file.readlines()
        global counter
        for i in lines:
            if i.find('mem=')!=-1:
                # add memory
                #if hasattr(atoms,"mem"):
                mod_lines.append('%mem={}MB\n'.format(self.memory))
                #else:
                #    mod_lines.append(i)
                continue
            if i.find('nprocshared=')!=-1:
                # add processors
                #if hasattr(atoms,"nproc"):
                mod_lines.append('%nprocshared={}\n'.format(self.nprocshared))
                #else:
                #    mod_lines.append(i)
                continue
            if i.find('PDBName')==-1:
                mod_lines.append(i)
                continue
            else:
                spl_line = re.findall(r'\S+|\s+', i)
                frozen.append(int(spl_line[3]))
                spl_line[3] = str(0)
                spl_line[5] = str(atoms.positions[c][0])
                spl_line[7] = str(atoms.positions[c][1])
                spl_line[9] = str(atoms.positions[c][2])
                c+=1
                mod_line = ''.join(spl_line)
                mod_lines.append(mod_line)
        open_file.close()
        directory_name = os.path.dirname(label)
        if directory_name:
            image_num = int(directory_name.split("-")[-1])
            counter[image_num]+=1
            if not os.path.exists(directory_name):
                os.makedirs(directory_name)
            label = label + '-{}'.format(int(counter[image_num]))
        else:
            counter[0]+=1
            label = label + '-{}'.format(int(counter[0]))

        self.prefix = os.path.basename(label)
        open_file = open(label + '.com', "wt")
        open_file.writelines(mod_lines)
        open_file.close()
        self.current_label = label
        self.frozen = frozen
                

    def read_results(self):
        output = read_gaussian_out(self.current_label + '.log')
        modify_forces(output, self.frozen)
        open_file = open(self.current_label + '-data.xyz', "wt")
        open_file.writelines(str(output.calc.results))
        open_file.close()
        self.calc = output.calc
        self.results = output.calc.results


