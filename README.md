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
 * In some systems zlib might be required. For example, in ubuntu one can get it by using `sudo apt install zlib1g-dev`. On a mac `brew install zlib` reportedly works.
 * `gcc`, `g++`, `gfortran` compilers
 * `cmake`
 
## Requirements [CYGWIN 2.900]
 * cmake
 * gcc-core
 * gcc-gfortran (gfortran)
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
 
It is now assumed that Ipopt from the coin-or project has been previously compiled with **HSL**.
This will enable **MA57** as the default linear solver.
However, it is crucial to verify that the coinor libraries have been generated in standard locations.
Namely, `/usr/local/lib`, and the libraries include: 

 * `coinasl`
 * `coinmetis`
 * `coinhsl`
 
Depending on your platform, these are typically named libcoinX.so or libcoinX.dll.a, where X is equal to asl, metis or hsl.

## Installation
 1. Visit the [Ipopt](https://github.com/coin-or/Ipopt) repository, follow the instructions and make sure it is compiled, with all the thirdparty libraries.
 2. At the root directory use cmake to generate the makefile e.g. `cmake CMakeLists.txt`
 3. Run `make`
 4. Check the bin directory to find the `k_aug` executable
 5. (*Windows*) add the LAPACK library directory to the PATH (*Cygwin* typically `/usr/lib/lapack`, as well the `libcoinX.dll.a` files (typically located at `/usr/local/lib/`)

## Installation (Mac OS X)
 1. Install IPOPT with the HSL libraries e.g. MA57. Please follow the instructions from the [documentation](https://coin-or.github.io/Ipopt/INSTALL.html)

Normally this will generate libraries in the `/usr/local/lib`, these are named `libcoinX.dylib` where X=asl,hsl, etc.

 2. Find out where libgfortran.dylib is located, tipically `/usr/local/Cellar/gcc/`... or `/usr/local/opt/gcc/lib/gcc/X/lib`.
Note that this will depend on the version of Mac OS X, gcc, etc.

 3. Put the location of gfortran in the line `80` of the `CMakeLists.txt` **after** `HINTS`.

 4. If you have libraries with different names other than `libcoinX.dylib` or different locations, make sure these are reflected in the lines `73-80` of the `CMakeLists.txt` file.
 5. Run `cmake .`, then `make` if successful, find your executables in the `bin` directory.
 6. Enjoy.

## Known issues
 * AMPL can not recognize command line options
 * Mac os is currently supported (I think!).
 
So far, dozens of times tested
[.](https://giphy.com/gifs/kSlJtVrqxDYKk/html5)

`k_aug` is an essential part of the NMPC-MHE framework(caprese). Written by David Thierry 2020, under BSD 3-Clause license.


