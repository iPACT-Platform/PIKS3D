### Short cuts (alias)

Put add this line on your ~/.bashrc

```bash
alias dvm3d='export PRJ=${HOME}/mpi_dvm3d; source ${PRJ}/conf/bashrc singleSphere'
```
NOTE: the singleSphere it the deafult run case.

Any time your login and needs to work on the project, type `dvm`. 
Then the aliases and related environments are set for you. 
For the compelete list of environment and aliases see `conf/bashrc`

### How to compile for the first time and the re-compilation

Supposing your are current in the project directory (the one with src and config direcories in it)

#### For the Release mode compilation

```bash
mkdir build
cd build
cmake ..
make -j 6
```

This will do the Release mode compilation. 
Anytime you make some modification on the source code, and you want to re-compile, just cd to the `build` directory, and make:

```
cd $BUILD
make -j 6
```

or simplely,

```bash
make -j 6 -C $BUILD
```

or even more simply using the alias we have set,

```bash
mb
```
which means make for build. The actural command mb issues can be seen by `alias mb`. 

To disable use of OpenMP, you can add the flag `-DDVM_OPENMP=OFF` to the
CMake command.

#### For the Debug mode compilation
For Debug build, similarly, start from the project directory

```bash
mkdir debug
cmake .. -DCMAKE_BUILD_TYPE=Debug
make
```
Anytime you make some modifications on the source code, and you want to re-compile for the Debug mode, just cd to the `debug` directory,
```
cd $DEBUG
make -j 6
```

or

```bash
make -j 6 -C $DEBUG
```

or using the alias

```bash
md
```

### How to run it

```bash
mkdir -p run/yourCase
cd run/yourCase
cp $CONF/para.in .
# do some adjusting in para.in
# then prepare the iamge data file in the current case direcory
# set the num of OMP threads via the environment
export OMP_NUM_THREADS=2
# run the release mode exe.:
mpirun -np 8 $BUILD/dvm3d.x 
# or run the debug mode exe.:
mpirun -np 8 $DEBUG/dmv3d.x
```

### About serial version code (v14) in the conf/

This file is `conf/linearised_BGK_D3Q64_v14.f90`.
This code is used as validation comparing purpose. 
For a serial run, always copy this file to the run case directory, and modify and compile then run.
If there is something need to change permnently, e.g., redefine the absicss/equlibrium, please change this file and commit also.
