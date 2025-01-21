# Code Descriptions

## 1. `gaussian_ONIOM_calculator.py` (present in [`modules`](https://github.com/2253shubham/Gaussian_ONIOM_plus_ASE_for_QMMM/tree/main/modules))
This is a derived class of the ASE Gaussian calculator that facilitates interfacing Gaussian ONIOM input and output files with ASE. 

### Key Features:
- Takes a Gaussian ONIOM input file template (examples in [`sample_input_files`](https://github.com/2253shubham/Gaussian_ONIOM_plus_ASE_for_QMMM/tree/main/sample_input_files)) and modifies atomic positions to generate new input files for the same structure.
- Assumes no modification is needed for parameters like functionals, basis set, symmetry options (specified in the link0 section).
- Ensures the QM and MM region distribution in the ONIOM input file remains fixed for the entire calculation.

**Note**: While creating the ONIOM input file, ensure all calculation parameters are correctly set, as these will not be modified during runtime.

---

## 2. `geo_opt.py`
This code performs geometry optimization of a given structure using QM/MM. It implements the `DampedDynamics` optimization class from ASE, though it can be modified.

### Inputs:
1. Input structure for optimization.
2. Reference Gaussian ONIOM template files.
3. Labels for Gaussian output files.
4. Tolerance for force convergence criteria.
5. Trajectory file name to store the converged trajectory (e.g., `.xyz` format).
6. Memory and processor count (`memory` and `nproc`).

### Output:
- The converged trajectory of the geometry optimization.

### Running on HPC Machines:
- Use the script `sub-opt-anvil.script` in [`bash_scripts`](https://github.com/2253shubham/Gaussian_ONIOM_plus_ASE_for_QMMM/tree/main/bash_scripts). Modify it as needed for other HPC machines.
- Setup:
  - Create a folder with:
    - Input structure file.
    - Linked codes (`gaussian_ONIOM_calculator.py`, `geo_opt.py`).
    - Scripts (`sub-opt-anvil.script`).
    - Reference Gaussian template files.
    - Optional forcefield parameter files (e.g., `oplsaa.prm`) if not using UFF.
    - A `ref_files.txt` in [`parameters`](https://github.com/2253shubham/Gaussian_ONIOM_plus_ASE_for_QMMM/tree/main/parameters) listing the Gaussian template file names.
  - Run the code via the command in `sub-opt-anvil.script` or submit it as a batch job.

**Note**: Allocate appropriate resources (`memory` and `nproc`) based on the Gaussian calculations.

---

## 3. `create_intermediates.py`
This code extracts intermediate structures from a converged trajectory, useful for creating input trajectories for NEB calculations.

### Inputs:
- Trajectory file for extraction.
- Reaction type.
- Number of intermediates to extract.
- Directory to store the intermediates.

---

## 4. `find_ts.py`
This code performs transition state (TS) searches using the NEB method.

### Inputs:
1. Reference Gaussian ONIOM files.
2. Labels for Gaussian output files.
3. Number of images in NEB calculations (excluding reactant and product states).
4. Number of guide images in the NEB trajectory (excluding reactant and product states).
5. Optimized reactant and product structure files.
6. Processors per image.
7. Memory and processor count (`memory` and `nproc`).

### Running on HPC Machines:
- Use `sub-neb.script` or `sub-neb-anvil.script` in `bash_scripts`. Modify the latter for other HPC machines.
- Setup:
  - Create a folder with:
    - Optimized reactant and product structure files.
    - Intermediate structure files generated using `create_intermediates.py`.
    - Linked codes (`gaussian_ONIOM_calculator.py`, `find_ts.py`).
    - Scripts (`sub-neb.script`, `sub-neb-anvil.script`).
    - Reference Gaussian ONIOM template files.
    - Optional forcefield parameter files (e.g., `oplsaa.prm`) if not using UFF.
    - `ref_files.txt` in the `parameters` folder listing Gaussian template file names.
    - `specs.txt` in the `parameters` folder specifying the number of images for NEB.

- Submit the job as a batch process on the HPC, ensuring necessary resources (`memory` and `nproc`) are allocated.

---
