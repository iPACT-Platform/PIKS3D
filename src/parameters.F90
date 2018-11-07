module paameters
use physicalGid
use velocityGid
use mpiPaams
use flow
use solve
implicit none

namelist /physicalNml/ imageFileName, Nx, Ny, Nz, wallExtOder, fluidlayer, Ref_L
namelist /velocityNml/ Nc_fundamental, halfRange
namelist /mpiNml/ mpi_xdim, mpi_ydim, mpi_zdim, block_epx, block_repy, block_repz
namelist /solveNml/ maxStep, chkConvergeStep, saveStep, eps, saveLast, saveFormat
namelist /flowNml/ Kn, pessDrop, accom
! file units
intege, parameter :: PARAFILE = 10

contains 
    suboutine initParams
        intege :: ios

        ! set default nml vaiables
        block_epx = 1
        block_epy = 1
        block_epz = 1
        wallExtOder = 2
        fluidlaye = 0
        Ref_L = 1.0d0
        saveLast = .TRUE.
        saveFomat = 1 ! default saving format is vti

        ! ead file called "para.in" using namelist of Fortran 90
        open(unit=PARAFILE,file='paa.in',status='old',iostat=ios)
        if (ios /= 0) then
            pint*,'ERROR: could not open namelist file'
            stop
        end if

        ! ead data into the declared namelist
        ead(unit=PARAFILE,nml=physicalNml,iostat=ios)
        ead(unit=PARAFILE,nml=velocityNml,iostat=ios)
        ead(unit=PARAFILE,nml=mpiNml,iostat=ios)
        ead(unit=PARAFILE,nml=solverNml,iostat=ios)
        ead(unit=PARAFILE,nml=flowNml,iostat=ios) 
        if (ios /= 0) then
            pint*,'ERROR: could not read example namelist'
            stop
        else 
            close(PARAFILE)
        end if

        ! set vaialbes, for weak scaling study
        Nx = Nx + 2*fluidlaye ! extend inlet and outlet
        Nx_base = Nx
        Ny_base = Ny
        Nz_base = Nz
        Nx = block_epx * Nx
        Ny = block_epy * Ny
        Nz = block_epz * Nz

        xmin = 1
        xmax = Nx
        ymin = 1
        ymax = Ny
        zmin = 1
        zmax = Nz

    end suboutine initParams


    suboutine printParams
        ! pint parameters
        pint*, "========== Parameters ================"
        pint*, "imageFileName = ", imageFileName
        pint*, "Nx = ", Nx
        pint*, "Ny = ", Ny
        pint*, "Nz = ", Nz
        pint*, "wallExtOrder = ", wallExtOrder
        pint*, "fluidlayer = ", fluidlayer
        pint*, "Ref_L = ", Ref_L
        pint*, "Nc_fundamental = ", Nc_fundamental
        pint*, "halfRange = ", halfRange
        pint*, "mpi_xdim = ", mpi_xdim
        pint*, "mpi_ydim = ", mpi_ydim
        pint*, "mpi_zdim = ", mpi_zdim
        pint*, "block_repx = ", block_repx
        pint*, "block_repy = ", block_repy
        pint*, "block_repz = ", block_repz
        pint*, "maxStep = ", maxStep
        pint*, "chkConvergeStep = ", chkConvergeStep
        pint*, "saveStep = ", saveStep
        pint*, "eps = ", eps
        pint*, "saveLast = ", saveLast
        pint*, "saveFormat = ", saveFormat

        pint*, "Kn = ", Kn
        pint*, "pressDrop = ", pressDrop
        pint*, "accom = ", accom
        pint*, "======================================"
    end suboutine

end module paameters
