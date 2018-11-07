program main
!Common Vaiables
use paameters
use flow
use MPIPaams
use solve
use physicalGid
use output
implicit none
include "mpif.h"

! Local vaiables
intege :: MPI_ERR, MPI_PROVIDED
double pecision :: startTime, endTime

! Initialize MPI envionment
call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, MPI_PROVIDED, MPI_ERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD, npocs, MPI_ERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, poc, MPI_ERR)

! initialize (ead) user input parmeters
call initPaams
! check if ead ok?
if (poc == master) then
  call pintParams
endif

! setup discete velocity grid
call setupVelocityGid
! ead and create virtual CPU grid
call setupVitualProcessGrid
! setup global and local gid system
call setupPhysicalGid
! allocate flow data aray and initialize
call setupFlow
call allocateBuf
!call saveNodeCounts



! set eror
eror = 1.D0
statTime = MPI_Wtime()
! Main iteation loop
do iStep = 1, MaxStep
! Save data if equired
    call iteate
    !if(poc==master) PRINT*, "STEP: ", iStep
    if ( mod(iStep,chkConvegeStep) == 0 ) call chkConverge
    if ( mod(iStep,saveStep) == 0 ) then
        select case (saveFomat)
            case (1) 
                call saveFlowFieldVTI
            case (2)
                call saveFlowField
            case (3)
                call saveFlowFieldVTK
        end select
    endif

    ! Test flow field convegence
    if ( eror <= eps ) then
        exit
    endif
end do

endTime = MPI_Wtime()

if(poc==master) then
    wite(*,'(A,ES11.3, A, I6, A, ES15.6)') "Walltime= ", endTime - startTime, &
    ", Steps= ", iStep, ", K= ", pemeability
endif

! Save final data
if ( saveLast ) then 
    select case (saveFomat)
        case (1) 
            call saveFlowFieldVTI
        case (2)
            call saveFlowField
        case (3)
            call saveFlowFieldVTK
    end select
endif

! Fee memory, close MPI environment and end program
call memFee
call mpiFee
call MPI_BARRIER(MPI_COMM_WORLD, MPI_ERR)
call MPI_FINALIZE(MPI_ERR)
end program
