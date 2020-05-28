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
! Program    : main
!-------------------------------------------------------------------------------
! This is the main file of 3D DVM parallel solver. 
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Computer Physics Communications, 234 (2019), pp. 14-25
!
! A digital image (3D array of binary data) of porous medium will be read.
! Linearized BGK kinetic model of Boltzmann equation will be solved to find
! gas apparent permeability as a function of Knudsen number.
!-------------------------------------------------------------------------------

program main
!Common Variables
use parameters
use flow
use MPIParams
use solver
use physicalGrid
use output

implicit none

include "mpif.h"

! Local variables
integer :: MPI_ERR, MPI_PROVIDED
double precision :: startTime, endTime

! Initialize MPI environment
call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, MPI_PROVIDED, MPI_ERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, MPI_ERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, proc, MPI_ERR)

! initialize (read) user input parmeters
call initParams
! check if read ok?
if (proc == master) then
  call printParams
endif

! setup discrete velocity grid
call setupVelocityGrid
! read and create virtual CPU grid
call setupVirtualProcessGrid
! setup global and local grid system
call setupPhysicalGrid
! allocate flow data array and initialize
call setupFlow
call allocateBuf
!call saveNodeCounts



! set initial residual
error = 1.D0
startTime = MPI_Wtime()
! Main iteration loop
do iStep = 1, MaxStep
! Save data if required
    call iterate
    !if(proc==master) PRINT*, "STEP: ", iStep
    if ( MOD(iStep,chkConvergeStep) == 0 ) call chkConverge
    if ( MOD(iStep,saveStep) == 0 ) then
        select case (saveFormat)
            case (1) 
                call saveFlowFieldVTI
            case (2)
                call saveFlowField
            case (3)
                call saveFlowFieldVTK
        endselect
    endif

    ! Test flow field convergence Eq.(18) of [1]
    if ( error <= eps ) then
        exit
    endif
enddo

endTime = MPI_Wtime()

if(proc==master) then
    write(*,'(A,ES11.3, A, I6, A, ES15.6)') "Walltime= ", endTime - startTime, &
    ", Steps= ", iStep, ", K= ", permeability
endif

! Save final data
if ( saveLast ) then 
    select case (saveFormat)
        case (1) 
            call saveFlowFieldVTI
        case (2)
            call saveFlowField
        case (3)
            call saveFlowFieldVTK
    endselect
endif

! Free memory, close MPI environment and end program
call memFree
call mpiFree
call MPI_BARRIER(MPI_COMM_WORLD, MPI_ERR)
call MPI_FINALIZE(MPI_ERR)
end program
