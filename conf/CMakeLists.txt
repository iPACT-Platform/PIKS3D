cmake_minimum_required(VERSION 3.2)

project(piks3d LANGUAGES Fortran)

option(DVM_OPENMP "Enable OpenMP" ON)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -fimplicit-none")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fbounds-check -ftracdback -Wall")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -implicitnone")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -check all -traceback")
else()
  MESSAGE(WARNING "Unknown Fortran compiler - default flags")
endif()

if(DVM_OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

find_package(MPI REQUIRED)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})


file(GLOB_RECURSE sources  src/*.F90)

add_executable(piks3d ${sources})