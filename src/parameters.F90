!/------------------------------------------------------------\
!                                                             |
!  $$$$$$$\ $$$$$$\ $$\   $$\  $$$$$$\   $$$$$$\  $$$$$$$\    |
!  $$  __$$\\_$$  _|$$ | $$  |$$  __$$\ $$ ___$$\ $$  __$$\   |
!  $$ |  $$ | $$ |  $$ |$$  / $$ /  \__|\_/   $$ |$$ |  $$ |  |
!  $$$$$$$  | $$ |  $$$$$  /  \$$$$$$\    $$$$$ / $$ |  $$ |  |
!  $$  ____/  $$ |  $$  $$<    \____$$\   \___$$\ $$ |  $$ |  |
!  $$ |       $$ |  $$ |\$$\  $$\   $$ |$$\   $$ |$$ |  $$ |  |
!  $$ |     $$$$$$\ $$ | \$$\ \$$$$$$  |\$$$$$$  |$$$$$$$  |  |
!  \__|     \______|\__|  \__| \______/  \______/ \_______/   |
!                                                             |
!       Parallel Image-based Kinetic Solver (3D)              |
!\------------------------------------------------------------/

! The MIT License (MIT)

! Copyright (c) 2017-2020 
!   Lianhua Zhu <zhulianhua121@gmail.com> 
!   and Minh-Tuan Ho <minhtuanho.vn@gmail.com>
!   and Yonghao Zhang <y.h.zhang168@gmail.com>

! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!-------------------------------------------------------------------------------
! module    : parameters
!-------------------------------------------------------------------------------
! This is a module reading simulation parameters specified by users
! in "para.in" file for 3D DVM solver.
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Computer Physics Communications, 234 (2019), pp. 14-25
!
!! @param imageFileName : file name of the digital porous image
!! @param Nx, Ny, Nz    : number of voxels in x, y, z direction
!! @param wallExtOrder  : order of extrapolation scheme for incomming f at solid 
!! @param fluidLayer 	: number of fluid layers will be added at (each) inlet, outlet
!!						  allows periodic BC available (see page 21 of Ref.[1])
!! @param Ref_L 		: reference length defined by Nx/Ny (assume that Ny=Nz)
!! @param Nc_fundamental: half of the number of the discrete velocities in one axis
!!						  The total number of discrete velocities is Nv=(2*Nc_fundamental)^3
!! @param halfRange 	: half-range or full-range Gauss-Hermit is chosen
!!						  halfRange=T means half-range, halfRange=F means full-range
!! @param mpi_xdim, mpi_ydim, mpi_zdim : number of MPI (processes) subdomains in x,y,z direction, Section 5.2 of Ref.[1]
!! @param maxStep 		: maximum number of iterations
!! @param chkConvergeStep : iteration interval for checking convergence criteria Eq.(18) of Ref.[1]	
!! @param saveStep 		: iteration interval for exporting output files	
!! @param eps 			: convergence tolerance Eq.(18) of Ref.[1]	
!! @param saveLast 		: whether data from the last iteration is save or not
!! @param saveFormat 	: flow field format (1) VTI (2) Tecplot (3) VTK	
!! @param Kn		 	: Knudsen number defined in Eq.(1) of [1]	
!! @param pressDrop		: pressure drop applied at inlet/outlet Eq.(10) of [1]
!! @param accom		 	: tangential momentum accommodation coefficient (TMAC)  Eq.(8) of [1]				
!-------------------------------------------------------------------------------

module parameters
use physicalGrid
use velocityGrid
use mpiParams
use flow
use solver
implicit none

namelist /physicalNml/ imageFileName, Nx, Ny, Nz, wallExtOrder, fluidlayer, Ref_L
namelist /velocityNml/ Nc_fundamental, halfRange
namelist /mpiNml/ mpi_xdim, mpi_ydim, mpi_zdim, block_repx, block_repy, block_repz
namelist /solverNml/ maxStep, chkConvergeStep, saveStep, eps, saveLast, saveFormat
namelist /flowNml/ Kn, pressDrop, accom
! file units
integer, parameter :: PARAFILE = 10

contains 
    subroutine initParams
        integer :: ios

        ! set default nml variables
        block_repx = 1
        block_repy = 1
        block_repz = 1
        wallExtOrder = 2
        fluidlayer = 0
        Ref_L = 1.0d0
        saveLast = .TRUE.
        saveFormat = 1 ! default saving format is vti

        ! read file called "para.in" using namelist of Fortran 90
        open(unit=PARAFILE,file='para.in',status='old',iostat=ios)
        if (ios /= 0) then
            print*,'ERROR: could not open namelist file'
            stop
        end if

        ! read data into the declared namelist
        read(UNIT=PARAFILE,NML=physicalNml,IOSTAT=ios)
        read(UNIT=PARAFILE,NML=velocityNml,IOSTAT=ios)
        read(UNIT=PARAFILE,NML=mpiNml,IOSTAT=ios)
        read(UNIT=PARAFILE,NML=solverNml,IOSTAT=ios)
        read(UNIT=PARAFILE,NML=flowNml,IOSTAT=ios) 
        if (ios /= 0) then
            print*,'ERROR: could not read example namelist'
            stop
        else 
            close(PARAFILE)
        end if

        ! set varialbes, for weak scaling study
        Nx = Nx + 2*fluidlayer ! extend inlet and outlet
        Nx_base = Nx
        Ny_base = Ny
        Nz_base = Nz
        Nx = block_repx * Nx
        Ny = block_repy * Ny
        Nz = block_repz * Nz

        xmin = 1
        xmax = Nx
        ymin = 1
        ymax = Ny
        zmin = 1
        zmax = Nz

    end subroutine initParams


    subroutine printParams
        ! print parameters
        print*, "========== Parameters ================"
        print*, "imageFileName = ", imageFileName
        print*, "Nx = ", Nx
        print*, "Ny = ", Ny
        print*, "Nz = ", Nz
        print*, "wallExtOrder = ", wallExtOrder
        print*, "fluidlayer = ", fluidlayer
        print*, "Ref_L = ", Ref_L
        print*, "Nc_fundamental = ", Nc_fundamental
        print*, "halfRange = ", halfRange
        print*, "mpi_xdim = ", mpi_xdim
        print*, "mpi_ydim = ", mpi_ydim
        print*, "mpi_zdim = ", mpi_zdim
        print*, "block_repx = ", block_repx
        print*, "block_repy = ", block_repy
        print*, "block_repz = ", block_repz
        print*, "maxStep = ", maxStep
        print*, "chkConvergeStep = ", chkConvergeStep
        print*, "saveStep = ", saveStep
        print*, "eps = ", eps
        print*, "saveLast = ", saveLast
        print*, "saveFormat = ", saveFormat

        print*, "Kn = ", Kn
        print*, "pressDrop = ", pressDrop
        print*, "accom = ", accom
        print*, "======================================"
    end subroutine

end module parameters
