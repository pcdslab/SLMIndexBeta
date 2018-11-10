# SLMIndexBeta
The implementation of SLM-Index v0.7beta integrated into Comet-MS v2016.03 codebase

Welcome to the SLMIndexBeta wiki!


# What do you need?
1. make
2. cmake
3. libdivsufsort
4. Windows: MinGW

## Install libdivsufsort
1. Get libdivsufsort from here: https://github.com/y-256/libdivsufsort
2. Follow the instructions to install the libdivsufsort.

### Linux
Assuming <lds> is the path to libdivsufsort home,
1. Navigate to <lds>/build/lib/CMakeFiles/divsufsort.dir
2. Execute the command: `ar rcs libdivsufsort.a`
3. Copy the libdivsufsort.a to <slmindexbeta>/CometSearch/libdivsufsort
4. Also copy the divsufsort.h, config.h and lfs.h from the installed directory to <slmindexbeta>/CometSearch/libdivsufsort

### Windows
You will need Eclipse C/C++ IDE and MinGW GCC/G++ to build a SLM-Index compatible libdivsufsort.
1. open Windows Powershell
2. Navigate to <lds> and execute `mkdir build`
3. Navigate to <lds> and execute `cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="../../" -G "Eclipse CDT4 - MinGW Makefiles" ..`
4. Open Eclipse, right click anywhere in the "Project Explorer" and "Import".
5. Put the path to <lds>/build and import the project.
6. Right click on the project and build the project.
7. Copy the built lib/libdivsufsort.dll.a and paste as <slmindexbeta>/CometSearch/libdivsufsort/Windows/libdivsufsort.a
8. Also copy the include/lfs.h, config.h and divsufsort.h to <slmindexbeta>/CometSearch/libdivsufsort/Windows 

# How to build
1. Open the Cygwin shell or Git Bash shell.
1. Navigate to SLM-Index home directory <slmindexbeta>
2. Execute the following command: `make -j<JOBS>`
3. Note that in above command, you can skip the `-j <JOBS>` if working with Linux.

# How to Run Sample
1. Unzip the <slmindexbeta>/samples/rat.tar.gz -> ratNew2cleavages.fasta
2. Traverse to <slmindexbeta>/output/bin directory and execute: `cp ../../samples/comet.params ./`
3. Run: `./comet.exe ../../samples/*.ms2`
3. The results (PIN and PEPXML) files will be placed in <slmindexbeta>/samples directory
