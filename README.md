<img src="/docs/kauglogo.png" alt="Drawing" width="350px"/>

# k\_aug *(kay-ogg)*

## Basic feature list:

 * Sensitivity matrices with K*dsdp = [d2Ldxdp^T dcdp^T]
 * Reduced Hessian  matrices
 * Hessian and Jacobian computation via ASL.
 * `dot_sens` computations (for sensitivity)
 
*k\_aug* is compatible with Pyomo and AMPL via ASL. The main functionality can be used through suffixes. 

## Requirements
 * [CMAKE](https://cmake.org/) Version 3.5 or higher.
 * In some systems zlib might be required. For example, in ubuntu one can get it by using `sudo apt install zlib1g-dev`. On a mac `brew install zlib` reportedly works.
 * gcc compilers including gfortran
 
## Requirements [CYGWIN]
 * cmake
 * gcc-core
 * gcc-gfortran(gfortran)
 * g++
 * git (Use the CYGWIN version. DO NOT use the windows version!)
 * make
 * wget
 * zlib-devel


## Dependencies
 * [MC19](http://www.hsl.rl.ac.uk/download/MC19/1.0.0/a/) (HSL)
 * [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) (Script that downloads and configures is included)
 * [SCOTCH](https://www.labri.fr/perso/pelegrin/scotch/) (Script that downloads and configures is included)
 * [OPENBLAS](https://www.openblas.net/) (Script that downloads and configures is included)
 * [MUMPS](http://mumps.enseeiht.fr/) (Script that downloads and configures is included)
 * [ASL](https://ampl.com/resources/hooking-your-solver-to-ampl/) (Script that downloads and configures is included) 
 * [MC30](http://www.hsl.rl.ac.uk/catalogue/mc30.html) (HSL) (NOW OPTIONAL)
 * [Pardiso](https://pardiso-project.org/) (NOT SUPPORTED ANYMORE)

## (NEW!) Windows installation video with Cygwin.
[Link](https://youtu.be/LrBRu8N-Isk)

## Installation
 1. Be sure to have a steady internet connection. Enter thirdparty directory and run the get.X or help.X scripts in the following order: 
    1. `ASL`
    2. `OpenBlas`
    3. `Metis`
    4. `Scotch`
    5. `Mumps`
 2. Download and place the .tar.gz file of mc19 into the thirdparty/hsl/mc19 directory and then run `help.mc19`
 3. At the root directory use cmake to generate the makefile e.g. `cmake CMakeLists.txt`
 4. Run `make`
 5. Check the bin directory to find the `k_aug` executable
 6. (Windows) add the OpenBLAS directory to the PATH

## Known issues
 * AMPL can not recognize command line options
 * MUMPS will often try to use multiple cores. It is preferable to turn off this functionallity, set `OMP_NUM_CORES=1` .
 * Mac os is currently not supported.
 
So far, hundreds of times tested.
`k_aug` is an essential part of the NMPC-MHE framework(caprese). Written by David Thierry 2017, under BSD 3-Clause license.


