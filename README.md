seekr2
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/seekr2/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/seekr2/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/seekr2/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/seekr2/branch/master)

```
##########################################      ########
   ###   ######  ######  ### ###  #####      ####      ####
  #####   ##  #   ##  #   ##  #    ##  #    ##           ###
  ##      ##      ##      ## #     ##  #     ##           ##
   ##     ####    ####    ###      ##  #                 ###
    ##    ##      ##      ###     #####                ####
  #####   ##  #   ##  #   ## #     ##  #            ###
   ###   ######  ######  ###  ##  ###  #         ## 
##########################################    ##
         #######   #  #      ###           ##            #
       ##       #########       ##       ##  ####    #####
         ###      #  #   #######        ##       #######

Simulation-Enabled Estimation of Kinetic Rates - Version 2
```

## Overview
Fast and versatile multiscale milestoning to compute molecular thermodynamics
and kinetics.

Prepare and run milestoning calculations in the OpenMM, NAMD, and/or Browndye2
simulation engines for the purposes of obtaining the kinetics and 
thermodynamics of molecular processes such as: ligand-receptor 
binding and unbinding, membrane permeability, internal molecular dynamics, 
and many other situations.

This is only a quickstart guide to get SEEKR2 up and running as soon as
possible. To see more detailed instructions and tutorials, please see the
docs/ subfolder and OTHER RESOURCES HERE.

## Install

### Dependencies
Many of the dependencies for SEEKR2 will be installed automatically, but
some must be installed separately

#### OpenMM

OpenMM is recommended for the MD stage of SEEKR2.

SEEKR2 needs the seekr2plugin in order to use OpenMM for MD simulations. 
INSTRUCTIONS FOR seekr2plugin INSTALLATION COMING SOON...

Alternatively, NAMD2 may be used for MD if desired. See the NAMD2 section
below for installation of NAMD2.

#### Browndye2

SEEKR2 needs Browndye2 if BD simulations will be run (necessary for k-on
calculations). Please see BROWNDYE2 LINK HERE for Browndye installation
instructions.

#### NAMD2

If OpenMM is not desirable or available for the MD simulations, SEEKR2 may 
use NAMD2 in order to run MD simulations. NAMD2 is frequently already 
available on shared computing resources like clusters and supercomputers.
Not all SEEKR2 options may be available using NAMD2.

If you wish to install NAMD2 yourself, please see NAMD2 LINK HERE for 
installation of NAMD.

### Install SEEKR2
To install SEEKR2, clone this repository and install the package:

```
git clone https://github.com/seekrcentral/seekr2.git
cd seekr2
python setup.py install
```

### Testing SEEKR2
To test SEEKR2, run the following command in the seekr2/ directory:

```
python setup.py test
```

Additional continuous integration tests may be run from the
seekr2/seekr2/continuous_integration/ directory.

## Run

A SEEKR2 calculation needs a "Model Input File" to run. Several examples may
be found in seekr2/seekr2/data. Execute the following commands within the 
seekr2/seekr2/ directory to run a sample calculation on the host-guest system:

(It is assumed that both Browndye2 and the seekr2plugin for OpenMM have
already been installed).

```
python prepare.py data/sample_input_mmvt_openmm.xml
python run.py any ~/test_mmvt_openmm/model.xml
python converge.py ~/test_mmvt_openmm/model.xml
python analyze.py ~/test_mmvt_openmm/model.xml
```

### Important Options

In general, SEEKR2 programs can be run with the '-h' argument to see all
available options.

MORE HERE

## Authors and Contributors

The following people have contributed directly to the coding and validation
efforts of SEEKR2 (listed an alphabetical order of last name). 
Thanks also to everyone who has helped or will help improve this project by 
providing feedback, bug reports, or other comments.

* Rommie Amaro (principal investigator)
* Ilker Deveci (contributor)
* Hilliary Frank (contributor)
* Anand Ojha (developer)
* Andy Stokely (developer)
* Natalie Timms (contributor)
* Lane Votapka (lead developer)
* Jeff Wagner (contributor)

### Copyright

Copyright (c) 2021, Lane Votapka


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
