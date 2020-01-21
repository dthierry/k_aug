<img src="/docs/kauglogo.png" alt="Drawing" width="350px"/>

# k\_aug *(kay-ogg)*

## Basic feature list:

 * Sensitivity matrices with K*dsdp = [d2Ldxdp^T dcdp^T]
 * Reduced Hessian  matrices
 * Hessian and Jacobian computation via ASL.
 * `dot_sens` computations (for sensitivity)
 
*k\_aug* is compatible with Pyomo and AMPL via ASL. The main functionality can be used through suffixes. 

## Requirements
 * [Ipopt](https://github.com/coin-or/Ipopt) Version 3.13.X, with its thirdparty libraries.
 * [CMAKE](https://cmake.org/) Version 3.5 or higher.
 * In some systems zlib might be required. For example, in ubuntu one can get it by using `sudo apt install zlib1g-dev`. On a mac `brew install zlib` reportedly works.
 * gcc compilers including g++, gfortran
 
## Requirements [CYGWIN 2.900]
 * cmake
 * gcc-core
 * gcc-gfortran(gfortran)
 * gcc-g++
 * git
 * make
 * wget
 * zlib-devel
 * patch
 * openblas
 * lapack-devel
 * pkg-config


## Dependencies (NEW!)
 * [MC19](http://www.hsl.rl.ac.uk/download/MC19/1.0.0/a/) (HSL)
It is now assumed that Ipopt from the coin-or project has been previously compiled with **HSL**.
This will enable MA57 as the default linear solver.
However, it is crucial to verify that the coinor libraries have been generated in standard locations.
Namely, `/usr/local/lib`, and the libraries include: 

 * `coinasl`
 * `coinmetis`
 * `coinhsl`
 
Depending on your platform, these are typically named libcoinX.so or libcoinX.dll.a, where X is equal to asl, metis or hsl.

## Installation
 1. Visit the [Ipopt](https://github.com/coin-or/Ipopt) repository, follow the instructions and make sure it is compiled, with all the thirdparty libraries.
 1. Download and place the .tar.gz file of mc19 into the thirdparty/hsl/mc19 directory and then run `help.mc19`
 3. At the root directory use cmake to generate the makefile e.g. `cmake CMakeLists.txt`
 4. Run `make`
 5. Check the bin directory to find the `k_aug` executable
 6. (Windows) add the LAPACK library directory to the PATH, as well the libcoinX.dll.a files (typically located at /usr/local/lib/)

## Known issues
 * AMPL can not recognize command line options
 * Mac os is currently supported (I think!).
 
So far, dozens of times tested.
![DOZENS!](https://giphy.com/gifs/kSlJtVrqxDYKk/html5)

`k_aug` is an essential part of the NMPC-MHE framework(caprese). Written by David Thierry 2020, under BSD 3-Clause license.


