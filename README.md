# CONFIGURATION & COMPILATION

To customize the build process enter the directory with the source code
and run

    cmake

Setting up an alternative compiler in a different folder (here a build-
folder within the source) is possible with

    FC=gfortran cmake ..

Compilation is done by

    make

To build the html documentation type

    make doc

The code is distributed under the GNU General Public License - see the
accompanying LICENSE file for more details. So feel free to experiment
with this.

Copyright (C) 2006-2018
Tobias Illenseer <tillense@astrophysik.uni-kiel.de>
Manuel Jung <mjung@astrophysik.uni-kiel.de>
Jannes Klee <jklee@astrophysik.uni-kiel.de>
