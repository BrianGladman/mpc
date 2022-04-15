======================
A. A Note On Licensing
======================

Where files in this distribution have been derived from files licensed
under Gnu GPL or LGPL license terms, their headers have been preserved 
in order to ensure that these terms will continue to be honoured.  

Other files in this distribution that have been created by me for use
in building MPIR and MPFR using Microsoft Visual Studio 2013 are 
provided under the terms of the LGPL version 2.1

Running the MPFR tests automatically uses Python, which must hence be 
installed if you want to run them.  

============================================
B. Compiling MPC with the Visual Studio 2013
============================================

The VC++ project files are intended for use with Visual Studio 
2013 Professional, but they can also be used with Visual C++ 
2013 Express to build win32 applications.  When the new Windows
7 SDK is released it will allow the Express product to build for
x64 as well. 

Building MPFR
-------------

These VC++ build projects are based on MPIR 2.0 and MPFR-3.0.0. It
is assumed that MPIR has already been built and that the directories
containing MPIR and MPFR are at the same level in the directory 
structure:

    mpir
        build.vc10
            dll     MPIR Dynamic Link Libraries 
            lib     MPIR Static Libraries
            ....
    mpfr
        build.vc10
            dll     MPFR Dynamic Link Libraries
            lib     MPFR Static Libraries
            ....
    mpc
        build.vc10
            dll
            lib
            mpc_dll
            dll_tests
            mpc_lib
            lib_tests
            
The root directory name of the MPIR version that is to be used in 
building MPFR should be 'mpir' with any version number such as in
'mpir-3.0' removed.
 
The MPC source distribution should be obtained and expanded into the
MPC root directory (e.g. mpc). After this the build project files 
should be added so that the build.vc10 sub-directory is in the
MPC root directory as shown above.  

The root directory names for mpir and mpfr are assumed to be 'mpir' 
and 'mpfr' (this makes it easier to use the latest version of MPIR 
and MPFR without having to update MPIR and MPFR library names and
locations when new versions are released).
        
The build project for MPC, mpc.sln, contains four build projects:

    mpc_dll         build MPC dll
    dll_tests       multiple tests for the MPC dll
    mpc_lib         build MPC static library
    lib_tests       multiple tests for the MPC static library

each of which supports four configurations:

    win32 or x64
    release or debug

In building the multiple test projects for both the static and DLL 
library, the add_test_lib project needs to be built before all the 
other test projects are built.
    
If you wish to use the Intel compiler, you need to convert the build files
by right clicking on the MPC top level Solution and then selecting the 
conversion option.

========================================
C. Using MPC with the Visual Studio 2013
========================================

When using the static MPC library in an application, the default Microsoft
run-time libraries (CRT) used are the static ones.  If you wish to use the
DLL run-time libraries, you will need to rebuild MPC with these libraries
set in the project property pages.

The DLL build of MPC uses the MPFR and MPIR DLLs so these need to be
available to any application built with the MPC DLL. To build an application
that uses these DLLs it is necessaery to define the preprocessor symbol
__GMP_LIBGMP_DLL when the application is built.

===================
D. Acknowledgements
===================

My thanks to:

1. The GMP team for their work on GMP and the MPFR team for their work 
   on MPFR and MPC
2. Patrick Pelissier, Vincent Lefèvre and Paul Zimmermann for helping
   to resolve VC++ issues in MPFR.
3. The MPIR team for their work on the MPIR fork of GMP.
4  Jeff Gilcrist for his help in testing, debugging and improving the
   readme.txt file giving the build instructions
 
       Brian Gladman, February 2013
