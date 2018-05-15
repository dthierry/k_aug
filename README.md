<img src="/docs/kauglogo.png" alt="Drawing" width="350px"/>

# k\_aug *(kay-ogg)*

## Basic useful feature list:

 * Reduced sensitivity
 * Hessian and Jacobian computation via ASL
 * K*dsdp = [d2Ldxdp^T dcdp^T]
 * Reduced Hessian
 * *dot* computations (for sensitivity)
 
*k\_aug* is compatible with Pyomo and AMPL via ASL. The main functionality can be used through suffixes. 

## Requirements
 * [CMAKE] Version 3.5 or higher.
 * [MC19]() (HSL)
 * [METIS]() (Script that downloads and configures is included)
 * [SCOTCH]() (Script that downloads and configures is included)
 * [OPENBLAS]() (Script that downloads and configures is included)
 * [MUMPS]() (Script that downloads and configures is included)
 * ASL 
 * [MC30](http://www.hsl.rl.ac.uk/catalogue/mc30.html) (HSL) (NOW OPTIONAL)
 * [Pardiso](https://pardiso-project.org/) (NOT SUPPORTED ANYMORE)

## Installation
 -1 Enter thirdparty directory and run the get.X or help.X scripts in the following order: 
  -1 `ASL`
  -2 `OpenBlas`
  -3 `Metis`
  -4 `Scotch`
  -5 `Mumps`
 -2 Download and place the .tar.gz file of mc19 into the thirdparty/hsl/mc19 directory and then run `help.mc19`
 -3 At the root directory use cmake to generate the makefile e.g. `cmake CMakeLists.txt`
 -4 Run `make`
 -5 Check the bin directory to find the `k_aug` executable
 
 ### Notes on the installation.
  * Sometimes libz might be required. In Unix-like environments use the package manager to find it and install.


## Known issues
 * AMPL can not recognize command line options
 * 
So far, hundreds of times tested.
*k\_aug* is part of the NMPC-MHE framework. Under BSD 3-Clause license.

by David M Thierry
