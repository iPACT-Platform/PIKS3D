# Intel Fotran setting
#MPIFC=mpiifort
#FC=ifort

# GNU Fotran setting
MPIFC=mpif90
FC=gfortran

# Release mode flags
#LFLAGS=  -g -fpe0 -fpe-all=0 -check -debug -traceback
#CFLAGS= -c -g -fpe-all=0 -fpe0 -check -debug -traceback

#LFLAGS=   -fopenmp 
#CFLAGS= -c   -fopenmp -O2

# Release mode flags
LFLAGS= -fopenmp
CFLAGS= -c -fopenmp -O2

# Intel Fortran setting

#Object files
OBJS   = velocityGrid.o physicalGrid.o mpiParams.o flow.o \
		 solver.o main.o gaussHermite.o parameters.o output.o
		 #string_utility_module.o

#Link into an excutable
piks3d: $(OBJS)
	$(MPIFC) $(LFLAGS) $^ -o $@

#Compile the modules
main.o : main.F90 flow.o mpiParams.o solver.o parameters.o output.o
	$(MPIFC) $(CFLAGS) $<
output.o : output.F90 flow.o mpiParams.o solver.o physicalGrid.o
	$(MPIFC) $(CFLAGS) $<
parameters.o : parameters.F90 velocityGrid.o physicalGrid.o \
	mpiParams.o solver.o flow.o# string_utility_module.o 
	$(MPIFC) $(CFLAGS) $<
#string_utility_module.o : string_utility_module.F90
#	$(MPIFC) $(CFLAGS) $<
velocityGrid.o : velocityGrid.F90 gaussHermite.o
	$(MPIFC) $(CFLAGS) $<
physicalGrid.o : physicalGrid.F90
	$(MPIFC) $(CFLAGS) $<
flow.o : flow.F90 velocityGrid.o physicalGrid.o mpiParams.o
	$(MPIFC) $(CFLAGS) $<
mpiParams.o : mpiParams.F90 physicalGrid.o
	$(MPIFC) $(CFLAGS) $<
solver.o : solver.F90 velocityGrid.o physicalGrid.o flow.o mpiParams.o
	$(MPIFC) $(CFLAGS) $<
gaussHermite.o : gaussHermite.F90
	$(MPIFC) $(CFLAGS) $<

# Build options
all: piks3d

clean:
	rm -f ./*.o ./*.mod ./piks3d
