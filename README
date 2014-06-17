SHORT - parallel unified multipole algorithm for Electromagnetics (Puma-EM). 
Aim is providing surface Method of Moments for Electromagnetics, enhanced 
by using the Multilevel Fast Multipole Method. Code is parallelized and runs on 
Desktops and clusters.


1. Introduction
---------------

This is Puma-EM, a Parallel Unified Multipole Algorithm for Electromagnetics.

The aim of Puma-EM is to solve surface integral equations that arise in 
Computational Electromagnetics, by using boundary elements methods, namely 
the Method of Moments. The method is enhanced by the use of the Multilevel 
Fast Multipole Method, which expedites the matrix-vector multiplication 
required by the iterative algorithm.

The term "Unified" refers to the fact that the Multipole Engine could be 
used for acoustics or mechanics. "Unified" should therefore be seen as a 
goal to reach. Hopefully the Open-Sourcing of the code is a step towards 
that goal.

Puma-EM is distributed under the terms of the GNU General Public License v3. 
See COPYING for more information.

For any questions/bugs/requirements, mail <vandenbosch.idesbald@gmail.com>.


2. Capabilities
---------------

This code is currently capable of computing the bistatic Radar Cross Section (RCS),
monostatic RCS and monostatic (Sythetic Aperture Radar) SAR of PEC targets. 
These targets can include junctions (i.e. plates that intersect volumes for example),
therefore the code can solve complex geometries.

Puma-EM has already solved problems containing more than 40 million of variables on
a "modest" cluster (2 HP proliant ML350 machines, quad core and 16 GB RAM each), but 
it is able to solve around 1-2 million on a Desktop/laptop having 2 GBytes of RAM.

Please note that Puma-EM has been developed and used only on Linux machines. 
Windows (R) is not supported (yet). For running Puma-EM on Windows machines,
you should consider installing a virtual linux machine through VMware or any
other virtualization solution.

Puma-EM should work on *NIX/*BSD systems, including MacOS (R) 10.x, although
it has not been tested yet.


3. Install Puma-EM
------------------

First, extract the *.tar.gz file that you downloaded in your home directory.

Second, cd into Puma-EM, then open the guide.pdf file. You can install the necessary
files locally or, if you have root access to the machine, you can install it for all
the users. The needed libraries and programs are:
- gmsh, an open-source mesh generator
- g++ and gfortran compilers (works with intel compilers too)
- blitz++, a C++ container library
- OpenMPI for the parallel code
- scipy and numpy, Python scientific libraries
- mpi4py, which is mpi for Python
- matplotlib, a Python plotting library.

There are automated scripts for installing Puma for some distributions. Due to a 
lack of time, not all of them are fully tested, but they should work. Type 
./install.sh in a terminal and follow the instructions. If unsure, install as per
the guide.


4. Ready to Go
--------------

You are now ready to rock'n'roll!! 

Choose a master node, open a terminal, and simply type:

  $ ./run.sh

See the doc/EXAMPLES file for other examples, then play around with the
simulation_parameters.py file.


5. Code documentation
---------------------

If you want to generate the documentation of the code itself, type 
(doxygen and LaTeX required):

  $ make documentation


