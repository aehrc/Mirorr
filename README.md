Mirorr README
=============

The main web site for Mirorr is here: http://aehrc.github.io/Mirorr/

This file describes how to compile and test the Mirorr source code. The
contents of is file are described below. The main web site is likely to contain more up-to-date information.


CONTENTS
--------

- DESCRIPTION
- CITATION
- LICENSE
- PRE-REQUISITES - SUMMARY
- PRE-REQUISITES - DETAILED INSTRUCTIONS
- BUILD INSTRUCTION
- TESTING THE BUILD
- FREQUENTLY ASKED QUESTIONS (FAQ)


DESCRIPTION
-----------

> MIRORR: Multimodal Image Registration using blOck-matching and Robust Regression

This package contains the source code of Mirorr, which implement both the Mirorr and SymMirorr robust rigid/affine image registration methods presented in:

David Rivest-Henault, Nicholas Dowson, Peter B. Greer, Jurgen Fripp, and Jason Dowling "Robust inverse-consistent affine CT-MR registration in MRI-assisted and MRI-alone prostate radiation therapy." Medical Image Analysis (In press), 2015.

The permanent citable DOI link to the original source code used in the production of this paper is here: [dx.doi.org](http://dx.doi.org/10.4225/08/55372DE407418). The GitHub version has been updated and is (semi-)actively developed. 

This software has been developed by a team of researchers from [CSIRO](http://www.csiro.au/)'s [The Australian E-Health Research Centre](http://aehrc.com/). See AUTHORS.txt for more details.

To use SymMirorr, run the program as follows:  `mirorr --reg-mode symmetric [other program arguments]`
To use Mirorr, run the program as follows:     `mirorr --reg-mode classic   [other program arguments]`


### Detailed description ###

The mirorr program implements a robust multimodal image registration method that is based on local correlation computed using a block-matching approach, as described in the paper indicated above. Mirror is, for the moment, only concerned with global registration, that is, using either a rigid or an affine transformation model. When using the default parameter set (--reg-mode symmetric), the registration method benefits from a half-way space definition to gain inverse-consistency high degree. That means that the order of input images on the command line has a notably reduced effect on the end result, simplifying analysis, and increasing robustness.

This algorithm has been extensively tested for CT-MR, MR-MR and MR-PET registrations.


CITATION
--------

If you use this program for scientific research, we would appreciate if you could cite the paper and/or the code mentioned in the Description section, above.


LICENSE
-------

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty information see LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the above copyright notice for more information.


PRE-REQUISITES - SUMMARY
------------------------

The following packages need to be installed to build Mirorr:
 - CMake
 - ITK, v4.X recommanded.
 - Boost, components: program_options filesystem system


PRE-REQUISITES - DETAILED INSTRUCTIONS
--------------------------------------

The Mirror application is built using the CMake system. In addition, it depends on several other libraries:

* The Insight Toolkit, with powerful pipelining, multithreading and
  medical image IO/manipulation routines, and large helpful user community.
    http://www.itk.org/ITK/resources/software.html
  Note: only versions of ITK 4 or higher are supported.

* Several Boost C++ libraries: including program_options, shared_ptr,
  timer, array, serialization, property_tree
    http://www.boost.org/users/download/

These libraries were selected both because they are useful and because
they have BSD/MIT style licenses.

All of the libraries will need to be downloaded and installed. In each
case, their source can be downloaded from the urls listed above and
compiled based on their respective instructions. 

However it is probably more convenient to use your operating system's
package manager for this purpose. Instructions for Ubuntu are as follows:

1. First add the neurodebian repository to your package manager so you
can conveniently get hold of the latest ITK, as described here:

http://neuro.debian.net/ 
e.g. for Ubuntu 12.04: 
  ```bash
  wget -O- http://neuro.debian.net/lists/precise.au.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
  sudo apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 2649A5A9
  sudo apt-get update
  ```

e.g. for Ubuntu 14.04: 
  ```bash
  wget -O- http://neuro.debian.net/lists/trusty.au.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
  sudo apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 2649A5A9 
  sudo apt-get update
  ```

2. Install the cmake build system and its gui
`sudo apt-get install cmake cmake-curses-gui`

3. Install ITK
`sudo apt-get install libinsighttoolkit4-dev`

4. Install all the boost libraries
`sudo apt-get install libboost-all-dev`


BUILD INSTRUCTION
-----------------

> Assuming you are in the root directory of the package:

```bash
mkdir build; 
cd build; 
cmake .. -DCMAKE_BUILD_TYPE=Release; 
make
```


TESTING THE BUILD
-----------------

```bash
cd ../test; 
./test.sh
```


FREQUENTLY ASKED QUESTIONS (FAQ)
--------------------------------

- **Mirorr is slow.** Did you forget to define the CMake variable CMAKE_BUILD_TYPE to Release (CMAKE_BUILD_TYPE=Release)? Typical speedup is about 8x. In addition, you can double check that the program is really using all available CPU cores (or ask for more cores on the cluster). It is possible to make the program use a specific number of cores with --nthreads.

- **Mirorr is still slow.** Sometimes, the two highest resolution settings of the image pyramid can be vastly redundant. You can try to remove one with the -b and -c options.

- **I want more speed.** If you have a GPU that supports OpenCL (most recent NVIDIA or ATI), you can try our GPU implementation. Compile the program with (USE_OPENCL=ON) and turn on GPU by using the --use-gpu-bm switch.


