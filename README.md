# Fosite

Fosite is a generic framework for the numerical solution of hyperbolic conservation
laws in generalized orthogonal coordinates. Its main purpose is the simulation of
compressible flows in accretion disks. The underlying numerical solution method
belongs to the family of unsplit conservative finite volume TVD schemes. The method
is 2nd order accurate in space and uses high order Runge-Kutta and multistep schemes
for time evolution. In addition to the pure advection code several source terms have
been implemented including viscous diffusion and gravitational acceleration.

Fosite is written with object-oriented patterns in Fortran 2003 and follows the
Structure of Arrays ([SoA](https://en.wikipedia.org/wiki/AOS_and_SOA)) layout,
operating on generic field datatypes. This allows for high performance on
modern architectures ([SIMD](https://en.wikipedia.org/wiki/SIMD)). It is parallelized
and vectorized. The software is thereby optimized for the [NEC SX-Aurora
TSUBASA Vector Engine](https://www.nec.com/en/global/solutions/hpc/sx/vector_engine.html).

Additional information and help:  
**Website**: www.astrophysik.uni-kiel.de/fosite/  
**IRC-Chat**: #fosite on freenode (http://webchat.freenode.net?channels=%23fosite)


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

Copyright (C) 2006-2018  
Tobias Illenseer <tillense@astrophysik.uni-kiel.de>  
Manuel Jung <mjung@astrophysik.uni-kiel.de>  
Jannes Klee <jklee@astrophysik.uni-kiel.de>  
