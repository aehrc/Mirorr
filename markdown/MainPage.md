Mirorr is software which has been developed at [The Australian E-Health Research Centre](http://aehrc.com/) for **3D rigid/affine** medical image registration that implements a **robust** and **inverse-consistent** algorithm. Mirror is a command line application that is relatively easy to use. Being robust, it is suitable for both **mono- and multi-modal** applications. Inverse-consistency means that a transform computed from image A to image B is **precisely the inverse** of the one computed from image B to image A. That means that it is not necessary to manage situations like the one illustrated below, simplifying your pipeline.

![Non inverse-consistent registration](images/fig1.png)

More formally, Mirorr implements a robust multimodal image registration method that is based on local correlations computed using a block-matching approach. Both rigid and affine transformation models are available. By default, Mirorr benefits from a half-way space definition and is inverse-consistent for most practical purposes. That means that the order of input images on the command line has a notably reduced effect on the end result, simplifying analysis, reducing opportunities for user error, and increasing robustness. This algorithm has been extensively tested for CT-MR, MR-MR and MR-PET registrations, using pelvis and head/brain datasets. An experimental GPU implementation is available. 

This software has been developed by [David Rivest-Hénault](https://github.com/drhenault), [Nick Dowson](http://baconmockup.com/600/800/), and a team of researchers from [CSIRO](http://www.csiro.au/)'s [The Australian E-Health Research Centre](http://aehrc.com/), see AUTHORS.txt for more details, and is made available under the terms of a liberal open source license, see LICENSE.txt.

This is the main web page for CSIRO Mirorr. 


CONTENTS
========

- [REFERENCE AND CITATION](http://aehrc.github.io/Mirorr/#reference-and-citation)
- [LICENSE](http://aehrc.github.io/Mirorr/#license)
- [BUILDING THE PROGRAM](http://aehrc.github.io/Mirorr/#building-the-program)
- [USER MANUAL AND EXAMPLES](http://aehrc.github.io/Mirorr/#user-manual-and-examples)
- [FREQUENTLY ASKED QUESTIONS (FAQ)](http://aehrc.github.io/Mirorr/#frequently-asked-questions-faq)
- [REPORTING AN ISSUE](reporting-an-issue)


REFERENCE AND CITATION
======================

Mirorr is the reference implementation for the _Mirror_ & _SymMirorr_ methods described in:

David Rivest-Hénault, Nicholas Dowson, Peter B. Greer, Jurgen Fripp, and Jason Dowling. **"Robust inverse-consistent affine CT-MR registration in MRI-assisted and MRI-alone prostate radiation therapy."** Medical Image Analysis (2015), DOI: [10.1016/j.media.2015.04.014](http://dx.doi.org/10.1016/j.media.2015.04.014). 

The permanent citable DOI link to the original source code used in the production of this paper is here: [dx.doi.org](http://dx.doi.org/10.4225/08/55372DE407418). The GitHub version has been updated and is (semi-)actively developed. 

If you use this program for scientific research or other related work, we would appreciate if you could cite the paper and/or the DOI code mentioned above.


LICENSE
=======

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty information see LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the above copyright notice for more information.


BUILDING THE PROGRAM
====================

The README.md (a simple text file) contains detailed build instructions. Briefly, this software uses a regular CMake build. The following packages need to be installed to build Mirorr:

 - [CMake](http://www.cmake.org/).
 - [ITK](http://www.itk.org/). Version v4.X recommended (tested with: 4.1.0 to 4.7.1).
 - [Boost](http://www.boost.org/). At the minimum you needs the following components: `program_options filesystem system timer`
 - If you want to use the GPU implementation, you need a working OpenCL development installation (generally provided by NVIDIA CUDA or the ATI equivalent).

With any modern operating system, there should be no need to build these packages since they are readily available through your favourite package manager (I'm partial to apt-get).


USER MANUAL AND EXAMPLES
========================

Important usage note
--------------------

The **Mirorr** program implements two methods described in the original paper: _Mirorr_ and _SymMirorr_. Generally, you will want to use _SymMirorr_, which is the default behaviour of **Mirorr** (we acknowledge that this can be a bit counter-intuitive).

If you really want to use the non inverse-consistent _Mirorr_, run the program as follows:
```bash
mirorr --reg-mode classic [other program arguments]
```
Otherwise, to use _SymMirorr_, run the program as follows:
```bash
mirorr --reg-mode symmetric [other program arguments]  # or simply:
mirorr [other program arguments]
```

Getting help
------------

Mirorr as a relatively extensive help page, and this is where you should start searching:
```bash
mirorr --help
```
If you can't find what you're looking for, feel free to contact the corresponding author indicated in the paper.


Very quick start
----------------

There are two minimal examples provided in the test directory. You can run them to validate your build, and read the .sh file to get some inspiration.

```bash
cd ../test; 
./test.sh
./test-gpu.sh  # Will crash if your build has not been compiled with OpenCL support.
```


FREQUENTLY ASKED QUESTIONS (FAQ)
================================

- **Mirorr is slow.** Did you forget to define the CMake variable CMAKE_BUILD_TYPE to Release (CMAKE_BUILD_TYPE=Release)? Typical speedup is about 8x. In addition, you can double check that the program is really using all available CPU cores (or ask for more cores on the cluster). It is possible to make the program use a specific number of cores with --nthreads.

- **Mirorr is still slow.** Sometimes, the two highest resolution settings of the image pyramid can be vastly redundant. You can try to remove one with the -b and -c options.

- **I want more speed.** If you have a GPU that supports OpenCL (most recent NVIDIA or ATI), you can try our GPU implementation. Compile the program with (USE_OPENCL=ON) and turn on GPU by using the --use-gpu-bm switch.


REPORTING AN ISSUE
==================

Although the developers are perfect, the source code might not be so. If you encounter an issue, please report it in the [GitHub Page](https://github.com/aehrc/Mirorr) (you need a free GitHub account for that). Alternatively, if you are C++ literate, your can ask for a pull request. I will try to respond quickly. Finally, if none of the above works well for you, feel free to contact the corresponding author indicated in the paper.

