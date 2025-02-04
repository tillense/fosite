# Fosite

Fosite is a generic framework for the numerical solution of hyperbolic conservation
laws in generalized orthogonal coordinates. Its main purpose is the simulation of
compressible flows in accretion disks. The underlying numerical solution method
belongs to the family of unsplit conservative finite volume TVD schemes. The method
is 2nd order accurate in space and uses high order Runge-Kutta schemes for time evolution.
In addition to the pure advection code several source terms have been implemented to
account for, e.g., shear stresses and gravitational forces.

Fosite is written in Fortran 2008 using an object-oriented design. It follows the
Structure of Arrays ([SoA](https://en.wikipedia.org/wiki/AOS_and_SOA)) programming paradigm
to allow for high data throughput on modern high performances architectures
([SIMD](https://en.wikipedia.org/wiki/SIMD)). The code has been vectorized and optimized for
amd64 / intel x86_64 architectures and the [NEC SX-Aurora TSUBASA Vector Engine](https://www.nec.com/en/global/solutions/hpc/sx/vector_engine.html).
It has been verified to compile with [gfortran](https://gcc.gnu.org/fortran/) v8.4 and above,
the current intel fortran compiler ifx v2025.0 from the [intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html),
AMDs flang Fortran compiler included in [AOCC](https://www.amd.com/en/developer/aocc.html) v5.0,
NVidias nvfortran compiler included in the [NVidia HPC SDK](https://developer.nvidia.com/hpc-sdk) v24.11
and the [LLVM flang compiler](https://github.com/llvm/llvm-project) v19.

Parallelization is implemented using domain decomposition with MPI communication. It has
been verified to run with MPICH, OpenMPI and NECs MPI libraries.

Additional information and help:  
**Website**: www.astrophysik.uni-kiel.de/fosite/  

## Compiling the Serial Version
To customize the build process enter the directory with the source code
and run

    cmake <path_to_source>

In order to get information about available options type

    cmake -L <path_to_source>

Setting up an alternative compiler in a different folder (here a build-
folder within the source) is possible with

    FC=gfortran cmake <path_to_source>

Compilation is done by

    make

## Compiling the Parallel Version
In order to build the parallel version type

    cmake -DPARALLEL=ON <path_to_source>

Compilation is done by

    make

Finally run

    mpirun -n X programname

where X is the number of cores that should be used. In the standard output
the partitioning can be checked at "MESH-----> MPI partitioning: X1:X2:X3" where
X1*X2*X3 should be equal to X.


## Compile on NEC-SX Aurora
On NEC-SX Aurora vector machines there is a toolchain provided. Run

    cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain-SX-Aurora.cmake <path_to_source>

## Legal

The code is distributed under the GNU General Public License - see the
accompanying LICENSE file for more details. So feel free to experiment
with this.

Copyright (C) 2006-2025

Tobias Illenseer <tillense@astrophysik.uni-kiel.de>  
Manuel Jung <mjung@astrophysik.uni-kiel.de>  
Jannes Klee <jklee@astrophysik.uni-kiel.de>  
