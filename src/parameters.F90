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
integer, parameter :: PARAfile = 10

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
        open(unit=PARAfile,file='para.in',status='old',iostat=ios)
        if (ios /= 0) then
            print*,'ERRor: could not open namelist file'
            stop
        end if

        ! read data into the declared namelist
        read(unit=PARAfile,nml=physicalNml,iostat=ios)
        read(unit=PARAfile,nml=velocityNml,iostat=ios)
        read(unit=PARAfile,nml=mpiNml,iostat=ios)
        read(unit=PARAfile,nml=solverNml,iostat=ios)
        read(unit=PARAfile,nml=flowNml,iostat=ios) 
        if (ios /= 0) then
            print*,'ERRor: could not read example namelist'
            stop
        else 
            close(PARAfile)
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
