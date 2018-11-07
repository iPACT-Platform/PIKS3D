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
call MPI_COMM_SIZE(MPI_COMM_WorLD, nprocs, MPI_ERR)
call MPI_COMM_RANK(MPI_COMM_WorLD, proc, MPI_ERR)

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



! set error
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
        end select
    endif

    ! Test flow field convergence
    if ( error <= eps ) then
        exit
    endif
end do

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
    end select
endif

! Free memory, close MPI environment and end program
call memFree
call mpiFree
call MPI_BARRIER(MPI_COMM_WorLD, MPI_ERR)
call MPI_FINALIZE(MPI_ERR)
end program
