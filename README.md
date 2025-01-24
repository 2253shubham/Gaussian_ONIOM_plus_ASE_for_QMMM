# QMMM Calculations implementing Gaussian ONIOM with ASE

<p class="center-content"> 
  <img src="https://github.com/2253shubham/Gaussian_ONIOM_plus_ASE_for_QMMM/blob/main/docs/Gaussian_plus_ASE.png" alt=""/>
</p>

This repository contains the scripts to perform QMMM calculations and implement them in geometry optimization and transition state searches. Gaussian ONIOM model is used to perform QMMM. The theory on Gaussian ONIOM can be seen [here](https://gaussian.com/oniom/). 

Currently Gaussian inbuilt transition state (TS) search methods are very finicky and not reliable where there are many weak vibrational modes associated either with torsional motions and/or framework atom vibrations present along with the vibration mode associated with the transition step (eg - bond breaking / formation). To circumvent this, we integrated the Gaussian calculations with ASE. ASE (Atomic Simulation Environment) supports many robust optimization and TS searches such as NEB (nudged elastic band). ASE provides Gaussian calculator object, which allows to execute gaussian calculations on a python interface.


## Documentation

[Documentation](https://github.com/2253shubham/Gaussian_ONIOM_plus_ASE_for_QMMM/blob/main/docs/Documentation.md)



## Authors

- [@shubham](https://github.com/2253shubham)


## Environment Variables

You may need to have an installed Gaussian software or if not, install Gaussian (version - g16) in your local system (or where you plan to ruin this code). Following which, the required environment variables can be set by the followiong lines:
```bash
module spider gaussian      # list installed Gaussian versions
module load gaussian        # load default version
```

Required python version  - 3.9 or higher \
Download and install [requirements.txt](https://github.com/2253shubham/Gaussian_ONIOM_plus_ASE_for_QMMM/blob/main/requirements.txt) to install python library dependencies.
```bash
pip install -r /path/to/requirements.txt
```
