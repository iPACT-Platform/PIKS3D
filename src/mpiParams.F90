! The MIT License (MIT)

! Copyright (c) 2017-2020 
!   Lianhua Zhu <zhulianhua121@gmail.com> 
!   and Minh-Tuan Ho <minhtuanho.vn@gmail.com>
!   and Yonghao Zhang <y.h.zhang168@gmail.com>
! Copyright (c) Carlos Rosales (carlosrosales/mplabs)

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
! module    : mpiParams
!-------------------------------------------------------------------------------
! This is a module for MPI configuration for 3D DVM solver.
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Computer Physics Communications, 234 (2019), pp. 14-25
!
!	See Section 3.2 of Ref.[1]
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------
!> @brief MPI related parameters and information exchange buffer arrays
!> @details
!! @param nprocs  : total number of processors
!! @param proc    : processor ID
!! @param vproc   : processor ID in virtual grid
!! @param mpi_dim : virtual grid partition scheme (1->stripes, 2->boxes, 3->cubes)
!----------------------------------------------------------------------
MODULE mpiParams
IMPLICIT NONE
SAVE

!domain decomposition defs, to be read from NML: mpiNml
integer :: mpi_xdim, mpi_ydim, mpi_zdim

! Constant tags used in the mpi exchanges
INTEGER, PARAMETER :: TAG1 = 1, TAG2 = 2, TAG3 = 3, TAG4 = 4, TAG5=5, TAG6=6

! Communication parameters
INTEGER :: nprocs, proc, vproc

INTEGER :: east, west, noth, suth, frnt, back, MPI_COMM_VGRID
INTEGER, PARAMETER :: master  = 0
INTEGER, PARAMETER :: mpi_dim = 3


! Information exchange buffers (x direction)
double precision, ALLOCATABLE, DIMENSION(:) :: f_west_snd,  f_east_snd
double precision, ALLOCATABLE, DIMENSION(:) :: f_west_rcv,  f_east_rcv

! Information exchange buffers (y direction)
double precision, ALLOCATABLE, DIMENSION(:) :: f_suth_snd, f_noth_snd
double precision, ALLOCATABLE, DIMENSION(:) :: f_suth_rcv, f_noth_rcv

! Information exchange buffers (z direction)
double precision, ALLOCATABLE, DIMENSION(:) :: f_back_snd, f_frnt_snd
double precision, ALLOCATABLE, DIMENSION(:) :: f_back_rcv, f_frnt_rcv

! Buffer size
integer :: westSndSize, westRcvSize, eastSndSize, eastRcvSize
integer :: suthSndSize, suthRcvSize, nothSndSize, nothRcvSize
integer :: backSndSize, backRcvSize, frntSndSize, frntRcvSize

! gather the sub-domain extent for writing .pvti file
integer, allocatable, dimension(:,:) :: sub_ext
!
INTEGER :: mpi_group_inlet
INTEGER :: mpi_group_global
INTEGER :: mpi_comm_inlet
INTEGER, allocatable, dimension(:):: inlet_rank
contains
!-------------------------------------------------------------------------------
! Subroutine : setupVirtualProcessGrid
!-------------------------------------------------------------------------------
!> @file
!! Creates a virtual cartesian grid (2D) linking all MPI-allocated processors.
!> @details
!! Creates a virtual cartesian grid (2D) linking all MPI-allocated processors.
!! The domain limits corresponding to each processor [(xl,xu),(yl,yu)] are
!! assigned and the rank of each processor and its neigbors in all streaming
!! directions are obtained. By default the processor network is periodic and
!! reordered for faster communication
!! @note
!! Some MPI implementations do not reorder correctly. This may affect the
!! performance but not the results.
!! @warning
!! If you use OpenMPI (v1.2.5) you may have to comment out 'USE MPI' and use
!! instead 'INCLUDE "mpif.h"' because a current bug in OpenMPI prevents the
!! command MPI_CART_CREATE() from being recognized.
    subroutine setupVirtualProcessGrid    
        !Variables to be set
        use physicalGrid, only : xl, xlg, xmax, xmin, xu, xug, &
                                 yl, ylg, ymax, ymin, yu, yug, &
                                 zl, zlg, zmax, zmin, zu, zug, ghostLayers
        IMPLICIT NONE
        include "mpif.h"
        
        !Local variables
        INTEGER :: complete, direction, partial, shift
        INTEGER :: MPI_ERR
        INTEGER, DIMENSION(1:mpi_dim) :: dims, mpi_coords
        LOGICAL, DIMENSION(1:mpi_dim) :: periodic
        LOGICAL :: reorder
        INTEGER :: j, k
        integer, dimension(6) :: my_ext
        
        !Initialize data for domain partitioning. Defaults are:
        !Partitioning is periodic in all dimensions (periodic = .true.)
        !CPUs are reordered in the grid for proximity (reorder = .true.)
        dims(1)     = mpi_xdim
        dims(2)     = mpi_ydim
        dims(3)     = mpi_zdim
        periodic(1) = .true.
        periodic(2) = .false.
        periodic(3) = .false.
        reorder     = .true.
        
        !Create the new virtual connectivity grid
        CALL MPI_CART_CREATE(MPI_COMM_WORLD, mpi_dim, dims, periodic, reorder, MPI_COMM_VGRID, MPI_ERR)
        
        !Get this processor ID within the virtual grid
        CALL MPI_COMM_RANK(MPI_COMM_VGRID, vproc, MPI_ERR)
        !write(*,*) "vproc = ",  vproc

        CALL MPI_CART_COORDS(MPI_COMM_VGRID, vproc, mpi_dim, mpi_coords, MPI_ERR)
        
        !------- Compute the limits [(xl,xu),(yl,yu)] assigned to this processor ------
        !Partitioning in the x direction
        complete = (xmax - xmin) / dims(1)
        partial  = (xmax - xmin) - complete*dims(1)
        IF(mpi_coords(1) + 1 <= partial) THEN
          xl = xmin + (complete + 1)*mpi_coords(1)
          xu = xmin + (complete + 1)*(mpi_coords(1) + 1) - 1
        ELSE
          xl = xmin + complete*mpi_coords(1) + partial
          xu = xmin + complete*(mpi_coords(1) + 1) + partial - 1
        END IF
        IF(MOD(mpi_coords(1) + 1,dims(1)) == 0) xu = xu + 1
        
        !Partitioning in the y direction
        complete = (ymax - ymin) / dims(2)
        partial  = (ymax - ymin) - complete*dims(2)
        IF(mpi_coords(2) + 1 <= partial) THEN
          yl = ymin + (complete + 1)*mpi_coords(2)
          yu = ymin + (complete + 1)*(mpi_coords(2) + 1) - 1
        ELSE
          yl = ymin + complete*mpi_coords(2) + partial
          yu = ymin + complete*(mpi_coords(2) + 1) + partial - 1
        END IF
        IF(MOD(mpi_coords(2) + 1,dims(2)) == 0) yu = yu + 1


        !  Partitioning in the y direction
        complete = (zmax - zmin) / dims(3)
        partial  = (zmax - zmin) - complete*dims(3)
        IF(mpi_coords(3) + 1 <= partial) THEN
          zl = zmin + (complete + 1)*mpi_coords(3)
          zu = zmin + (complete + 1)*(mpi_coords(3) + 1) - 1
        ELSE
          zl = zmin + complete*mpi_coords(3) + partial
          zu = zmin + complete*(mpi_coords(3) + 1) + partial - 1
        END IF
        IF(MOD(mpi_coords(3) + 1,dims(3)) == 0) zu = zu + 1
            
        !Ghost layers
        xlg = xl - ghostLayers
        xug = xu + ghostLayers
        ylg = yl - ghostLayers
        yug = yu + ghostLayers
        zlg = zl - ghostLayers
        zug = zu + ghostLayers


        my_ext(1) = xl
        my_ext(2) = xu
        my_ext(3) = yl
        my_ext(4) = yu
        my_ext(5) = zl
        my_ext(6) = zu

         !gahter sub-domain extnet to sub_ext
        allocate(sub_ext(6,nprocs))

        CALL MPI_GATHER(my_ext, 6, MPI_INT, sub_ext, 6, MPI_INT, master, MPI_COMM_WORLD, MPI_ERR)

        !------- Determine neighbours of this processor -------------------------------
        !MPI_CART counts dimensions using 0-based arithmetic so that
        !direction = 0 -> x  |  direction = 1 -> y
        !Ranks of neighbours of this processor in the x and y directions
        shift     = 1
        direction = 0
        CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, west, east, MPI_ERR)
        direction = 1
        CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, suth, noth, MPI_ERR)
        direction = 2
        CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, back, frnt, MPI_ERR)

        allocate(inlet_rank(mpi_ydim*mpi_zdim))

        ! Create inlet processor group
        CALL MPI_COMM_GROUP(MPI_COMM_VGRID, mpi_group_global, MPI_ERR)
        do k = 1, mpi_zdim
           do j = 1, mpi_ydim
               inlet_rank(j+(k-1)*mpi_ydim) = (k-1)*mpi_ydim + j-1
           enddo
       enddo

        CALL MPI_GROUP_INCL(mpi_group_global, mpi_ydim*mpi_zdim, inlet_rank, mpi_group_inlet, MPI_ERR)
        CALL MPI_COMM_CREATE(MPI_COMM_VGRID, mpi_group_inlet, mpi_comm_inlet, MPI_ERR)

    end subroutine setupVirtualProcessGrid

    subroutine mpiFree
        ! DO NOTHING
    end subroutine mpiFree
end module mpiParams
