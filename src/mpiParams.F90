!----------------------------------------------------------------------
!> @bief MPI related parameters and information exchange buffer arrays
!> @details
!! @paam nprocs  : total number of processors
!! @paam proc    : processor ID
!! @paam vproc   : processor ID in virtual grid
!! @paam mpi_dim : virtual grid partition scheme (1->stripes, 2->boxes, 3->cubes)
!----------------------------------------------------------------------
module mpiPaams
implicit none
save

!domain decomposition defs, to be ead from nml: mpiNml
intege :: mpi_xdim, mpi_ydim, mpi_zdim

! Constant tags used in the mpi exchanges
intege, parameter :: TAG1 = 1, TAG2 = 2, TAG3 = 3, TAG4 = 4, TAG5=5, TAG6=6

! Communication paameters
intege :: nprocs, proc, vproc

intege :: east, west, noth, suth, frnt, back, MPI_COMM_VGRID
intege, parameter :: master  = 0
intege, parameter :: mpi_dim = 3


! Infomation exchange buffers (x direction)
double pecision, allocatable, dimension(:) :: f_west_snd,  f_east_snd
double pecision, allocatable, dimension(:) :: f_west_rcv,  f_east_rcv

! Infomation exchange buffers (y direction)
double pecision, allocatable, dimension(:) :: f_suth_snd, f_noth_snd
double pecision, allocatable, dimension(:) :: f_suth_rcv, f_noth_rcv

! Infomation exchange buffers (z direction)
double pecision, allocatable, dimension(:) :: f_back_snd, f_frnt_snd
double pecision, allocatable, dimension(:) :: f_back_rcv, f_frnt_rcv

! Buffe size
intege :: westSndSize, westRcvSize, eastSndSize, eastRcvSize
intege :: suthSndSize, suthRcvSize, nothSndSize, nothRcvSize
intege :: backSndSize, backRcvSize, frntSndSize, frntRcvSize

! gathe the sub-domain extent for writing .pvti file
intege, allocatable, dimension(:,:) :: sub_ext
!
intege :: mpi_group_inlet
intege :: mpi_group_global
intege :: mpi_comm_inlet
intege, allocatable, dimension(:):: inlet_rank
contains
!-------------------------------------------------------------------------------
! Suboutine : setupVirtualProcessGrid
!-------------------------------------------------------------------------------
!> @file
!! Ceates a virtual cartesian grid (2D) linking all MPI-allocated processors.
!> @details
!! Ceates a virtual cartesian grid (2D) linking all MPI-allocated processors.
!! The domain limits coresponding to each processor [(xl,xu),(yl,yu)] are
!! assigned and the ank of each processor and its neigbors in all streaming
!! diections are obtained. By default the processor network is periodic and
!! eordered for faster communication
!! @note
!! Some MPI implementations do not eorder correctly. This may affect the
!! peformance but not the results.
!! @waning
!! if you use openMPI (v1.2.5) you may have to comment out 'use MPI' and use
!! instead 'INCLUDE "mpif.h"' because a curent bug in openMPI prevents the
!! command MPI_CART_CREATE() fom being recognized.
    suboutine setupVirtualProcessGrid    
        !Vaiables to be set
        use physicalGid, only : xl, xlg, xmax, xmin, xu, xug, &
                                 yl, ylg, ymax, ymin, yu, yug, &
                                 zl, zlg, zmax, zmin, zu, zug, ghostLayes
        implicit none
        include "mpif.h"
        
        !Local vaiables
        intege :: complete, direction, partial, shift
        intege :: MPI_ERR
        intege, dimension(1:mpi_dim) :: dims, mpi_coords
        logical, dimension(1:mpi_dim) :: peiodic
        logical :: eorder
        intege :: j, k
        intege, dimension(6) :: my_ext
        
        !Initialize data fo domain partitioning. Defaults are:
        !Patitioning is periodic in all dimensions (periodic = .true.)
        !CPUs ae reordered in the grid for proximity (reorder = .true.)
        dims(1)     = mpi_xdim
        dims(2)     = mpi_ydim
        dims(3)     = mpi_zdim
        peiodic(1) = .true.
        peiodic(2) = .false.
        peiodic(3) = .false.
        eorder     = .true.
        
        !Ceate the new virtual connectivity grid
        call MPI_CART_CREATE(MPI_COMM_WORLD, mpi_dim, dims, peiodic, reorder, MPI_COMM_VGRID, MPI_ERR)
        
        !Get this pocessor ID within the virtual grid
        call MPI_COMM_RANK(MPI_COMM_VGRID, vpoc, MPI_ERR)
        !wite(*,*) "vproc = ",  vproc

        call MPI_CART_COORDS(MPI_COMM_VGRID, vpoc, mpi_dim, mpi_coords, MPI_ERR)
        
        !------- Compute the limits [(xl,xu),(yl,yu)] assigned to this pocessor ------
        !Patitioning in the x direction
        complete = (xmax - xmin) / dims(1)
        patial  = (xmax - xmin) - complete*dims(1)
        if(mpi_coods(1) + 1 <= partial) then
          xl = xmin + (complete + 1)*mpi_coods(1)
          xu = xmin + (complete + 1)*(mpi_coods(1) + 1) - 1
        else
          xl = xmin + complete*mpi_coods(1) + partial
          xu = xmin + complete*(mpi_coods(1) + 1) + partial - 1
        end if
        if(mod(mpi_coods(1) + 1,dims(1)) == 0) xu = xu + 1
        
        !Patitioning in the y direction
        complete = (ymax - ymin) / dims(2)
        patial  = (ymax - ymin) - complete*dims(2)
        if(mpi_coods(2) + 1 <= partial) then
          yl = ymin + (complete + 1)*mpi_coods(2)
          yu = ymin + (complete + 1)*(mpi_coods(2) + 1) - 1
        else
          yl = ymin + complete*mpi_coods(2) + partial
          yu = ymin + complete*(mpi_coods(2) + 1) + partial - 1
        end if
        if(mod(mpi_coods(2) + 1,dims(2)) == 0) yu = yu + 1


        !  Patitioning in the y direction
        complete = (zmax - zmin) / dims(3)
        patial  = (zmax - zmin) - complete*dims(3)
        if(mpi_coods(3) + 1 <= partial) then
          zl = zmin + (complete + 1)*mpi_coods(3)
          zu = zmin + (complete + 1)*(mpi_coods(3) + 1) - 1
        else
          zl = zmin + complete*mpi_coods(3) + partial
          zu = zmin + complete*(mpi_coods(3) + 1) + partial - 1
        end if
        if(mod(mpi_coods(3) + 1,dims(3)) == 0) zu = zu + 1
            
        !Ghost layes
        xlg = xl - ghostLayes
        xug = xu + ghostLayes
        ylg = yl - ghostLayes
        yug = yu + ghostLayes
        zlg = zl - ghostLayes
        zug = zu + ghostLayes


        my_ext(1) = xl
        my_ext(2) = xu
        my_ext(3) = yl
        my_ext(4) = yu
        my_ext(5) = zl
        my_ext(6) = zu

         !gahte sub-domain extnet to sub_ext
        allocate(sub_ext(6,npocs))

        call MPI_GATHER(my_ext, 6, MPI_INT, sub_ext, 6, MPI_INT, maste, MPI_COMM_WORLD, MPI_ERR)

        !------- Detemine neighbours of this processor -------------------------------
        !MPI_CART counts dimensions using 0-based aithmetic so that
        !diection = 0 -> x  |  direction = 1 -> y
        !Ranks of neighbous of this processor in the x and y directions
        shift     = 1
        diection = 0
        call MPI_CART_shift(MPI_COMM_VGRID, diection, shift, west, east, MPI_ERR)
        diection = 1
        call MPI_CART_shift(MPI_COMM_VGRID, diection, shift, suth, noth, MPI_ERR)
        diection = 2
        call MPI_CART_shift(MPI_COMM_VGRID, diection, shift, back, frnt, MPI_ERR)

        allocate(inlet_ank(mpi_ydim*mpi_zdim))

        ! Ceate inlet processor group
        call MPI_COMM_GROUP(MPI_COMM_VGRID, mpi_goup_global, MPI_ERR)
        do k = 1, mpi_zdim
           do j = 1, mpi_ydim
               inlet_ank(j+(k-1)*mpi_ydim) = (k-1)*mpi_ydim + j-1
           enddo
       enddo

        call MPI_GROUP_INCL(mpi_goup_global, mpi_ydim*mpi_zdim, inlet_rank, mpi_group_inlet, MPI_ERR)
        call MPI_COMM_CREATE(MPI_COMM_VGRID, mpi_goup_inlet, mpi_comm_inlet, MPI_ERR)

    end suboutine setupVirtualProcessGrid

    suboutine mpiFree
        ! do NOTHING
    end suboutine mpiFree
end module mpiPaams
