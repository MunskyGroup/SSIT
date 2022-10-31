# FSPCVODE_MEX

Codes for building Matlab interface with Sundials's CVODES library using Mex. These codes target only initial value problems of the form dx/dt = A(t)*x.

Prerequisites:
- CMake.
- Sundials.
- Matlab.

To build (using command line):
- Create an empty subfolder within the FSPCVODE_MEX folder called "build". 
- Change directory (via, e.g, cd) to "build".
- Type "cmake ..".
- Type "make". The C source files are compiled into shared libraries as well as MEX files.
- Make sure the MEX files are added to MATLAB's search path.
- The Mex solver is ready to use. For more information, see the comments in "FspCVodeMex.c".

Advanced options:
1) Build with OpenMP support:
- The library can optionally use Sundials' OpenMP support (provided that this was enabled for your Sundials build). To enable OpenMP, instead of "cmake..", type

cmake -DUSE_OPENMP=ON ..

Note that the only compiler supported by Matlab on MacOS does not support OpenMP. If you are building this library in MacOS, the OpenMP support will always be turned off.
