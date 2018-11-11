# SLMIndexBeta
The implementation of SLM-Index v0.7beta integrated into Comet-MS v2016.03 codebase

# What do you need?
1. make
2. cmake
3. libdivsufsort
4. Windows: MinGW

## Install libdivsufsort
1. Get libdivsufsort from here: https://github.com/y-256/libdivsufsort
2. libdivsufsort can be configured to use OpenMP by modifying CMakeLists.txt. 
3. Follow the instructions to install the libdivsufsort.

### Linux
Assuming /lds is the path to libdivsufsort home,
1. Navigate to /lds/build/lib/CMakeFiles/divsufsort.dir
2. Execute the command: `ar rcs libdivsufsort.a`
3. Copy the libdivsufsort.a to /slmindexbeta/CometSearch/libdivsufsort
4. Also copy the divsufsort.h, config.h and lfs.h from the installed directory to /slmindexbeta/CometSearch/libdivsufsort

### Windows
You will need Eclipse C/C++ IDE and MinGW GCC/G++ to build a SLM-Index compatible libdivsufsort.
1. open Windows Powershell
2. Navigate to <lds> and execute `mkdir build`
3. Navigate to <lds> and execute `cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="../../" -G "Eclipse CDT4 - MinGW Makefiles" ..`
4. Open Eclipse, right click anywhere in the "Project Explorer" and "Import".
5. Put the path to <lds>/build and import the project.
6. Right click on the project and build the project.
7. Copy the built lib/libdivsufsort.dll.a and paste as /slmindexbeta/CometSearch/libdivsufsort/Windows/libdivsufsort.a
8. Also copy the include/lfs.h, config.h and divsufsort.h to <slmindexbeta>/CometSearch/libdivsufsort/Windows 

# How to build SLM-Index Beta
1. Open the Cygwin shell or Git Bash shell.
1. Navigate to SLM-Index home directory <slmindexbeta>
2. Execute the following command: `make -j<JOBS>`
3. Note that in above command, you can skip the `-j <JOBS>` if working with Linux.

# How to Run samples
1. Unzip the /slmindexbeta/samples/rat.tar.gz -> ratNew2cleavages.fasta
2. Traverse to /slmindexbeta>/output/bin directory and execute: `cp ../../samples/comet.params ./`
3. Run: `./comet.exe ../../samples/*.ms2`
4. The results (PIN and PEPXML) files will be placed in <slmindexbeta>/samples directory

# How to search MS/MS data against SLM-Index
1. Digest the proteome database using Protein Digestion Simulator or OpenMS.
2. Remove the redundant peptide sequences in the digested database using DBToolkit.
3. Optional: The decoy database can be generated and appended to the target database using DBToolkit

## Edit the comet.params file
1. Add the path to the database file.
2. Configure the number and types of modifications that should be added to SLM-Index Beta.
3. Configure the digestion parameters for SLM-Index to validate the sequences in the database file.
4. Configure the result and output format.
4. Configure the capacity of modified peptides that could be added to SLM-Index Beta. 

### Note: 
1. If capacity is smaller than the number of modified peptides generated from the normal peptides, SLM-Index Beta will exit and throw a error message for the user to either increase the capacity or reduce the number/types of modifications to be added.
2. The max fragment charge: 3
3. Max digestion mass: 65kDa
4. The num_threads paramter controls the number of threads that can be used for indexing. The querying is restricted to 1 core.

# Please cite us if you use our work
For queries about SLM-Index, please contact: fsaeed@fiu.edu or mhaseeb@fiu.edu. Thanks
