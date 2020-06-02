


> Written with This tutorial simulates the rarefied gas flow in a fractal model, which described by the level-3 Sierpinski carpet.

<p align="center"> <a href="https://ibb.co/p2kG22Z"><img src="https://i.ibb.co/0qz8qqZ/sandPack.png" alt="sandPack" border="0" width="360"> </a> </p>

The black and white regions represent the fluid and solid parts respectively.



### 1. Compile the source codes

Suppose that your are currently in the project directory (the one with `src/`, `conf/`, `references/`, `tutorials/` directories in it) and a MPI environment is ready

```bash
cd src/
make
```
This will create executive file `piks3d` inside the `src/`  directory. This compilation step will not have to be repeated, except you modify the source codes. 

### 2. Run simulations

In the `tutorials` directory, there are a few subfolders for some example cases. In each subfolder of `tutorials/` folder, there are an image file and a `para.in`  file. 

Change your working directory to the current tutorial from `src/`
```bash
cd ../tutorials/sandPack/
```
To adjust simulation setup, you can modify the values of parameters in the `para.in` file and/or replace the image file.

Next, run the solver using the following command from your current working directory

```bash
OMP_NUM_THREADS=8 mpirun -np 4 ../../src/piks3d

```
which will use, for example, 4 MPI processes and 8 OpenMP threads for each MPI process.

### 3. Simulation outputs

All simulation outputs are stored in the working directory `tutorials/sandPack/`:
*  A `Results.dat` file logs the screen output, which contains the permeability values for all of the cases, and their convergence history. Below is and example of `Results.dat` file when the all simulations finish. From the left to right, columns for Knudsen number, dimensionless mass flow rate, dimensionless permeability, residual and iteration number are represented. 
```
   .....
```
* Flow field files, e.g. `Field.pvti`. 

### 4. How to adjust simulation parameters

All simulation parameters are set in the `para.in` file.  Below is an example of `para.in` file in the directory `tutorials/sandPack/`
```
&physicalNml
imageFileName = 'SandPack300^3.raw'
Nx = 300
Ny = 300
Nz = 300
wallExtOrder = 1
fluidLayer = 4
Ref_L = 1.d0
/

&velocityNml
Nc_fundamental = 4
halfRange = T
/

&mpiNml
mpi_xdim = 2
mpi_ydim = 2
mpi_zdim = 1
/

&solverNml
maxStep = 100
chkConvergeStep = 10
saveStep = 50
eps = 1.0d-7
saveLast = T
saveFormat = 1
/

&flowNml
Kn = 1.0d-1
pressDrop = 1.0d-1
accom = 1.d0
/
```
where
* Line starting with & means a parameter group, do not change!
* `imageFileName` is the file name of the image file, i.e. geometry flag (0 for fluid node, and 1 for solid node) file in binary format. Note that the format of the 3D raw image file is well illustrated by the standalone Fortran source code `conf/createSphere.F90`. See its header for detailed explanation.
* `Nx`, `Ny` and `Nz` is the number of nodes in x-, y- and z-direction in the spatial space, should be consistent to the image file.
* `wallExtOrder` is  the accuracy order of extrapolation at the wall nodes. can be [1|2|3]. 1 means first order; 2 means second order; 3 means use first order when the wall is one fluid node off the communication boundary, and use second order otherwise. 
* `fluidLayer` is the number of fluid-node layers padded at the inlet and outlet.
* `Ref_L` is the reference length, which is equal to `Nx/Ny` (suppose that the cross-section perpendicular to the streamline (x-)direction is square `Ny=Nz`).
* `Nc_fundamental` is number of discrete point in half of an axis.  Here 4 means a total of $(2x4)^3$ for 3D problems. Available values are  2, 4, 6, 8, 12, 16.
* `halfRange` can be [T|F]: choose half-range Gauss-Hermite or not, either T or F.
* `mpi_xdim`, `mpi_ydim` and `mpi_zdim` : these are the numbers of subdomains in x-, y- and z-direction in a frame of Cartesian even domain decomposition. The number of MPI processes corresponds to the multiple of `mpi_xdim` and `mpi_ydim` and `mpi_zdim`.  
* `maxSteps` is the maximum number of iterations.
* `chkConvergeStep` the frequency to check the convergence of the permeability.
* `saveStep` the frequency to save the flow fields.
* `saveFormat` the format to save data, can be [1|2|3]. 1 means in .pvti format, 2 means in Tecplot format, 3 means in vtk format.
* `eps` is value for the convergence criterion. The simulation will stop either the `eps` or the `maxSteps` is achieved, whichever comes first. Therefore, to achieve a reasonable `eps`, the maximum number of iterations `maxSteps` normally should be much larger than 100 as in this example.
* `saveLast` do we save the flow field after finishing running? Can be [T|F].
* `Kn` is Knudsen number.
* `pressDrop` is the non-dimensional pressure drop between the inlet and outlet. With the current linearized formulation and nondimensionalization, the dimensionless permeability are independent of positive`pressDrop`.  
* `accom` is the accommodation coefficient in the diffuse-specular model for the gas-solid interaction.

**Note**: 
More details on the simulation parameters can be found in Reference articles in the `references/` directory. 


