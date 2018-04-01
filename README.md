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

 * [MC30](http://www.hsl.rl.ac.uk/catalogue/mc30.html) (HSL)
 * [Pardiso](https://pardiso-project.org/)
 * ASL
 * Lapack
 * Blas

## Known issues
 * AMPL cannon recognize command line options
 * 
So far, hundreds of times tested.
*k\_aug* is part of the NMPC-MHE framework. Under BSD 3-Clause license.

by David M Thierry
