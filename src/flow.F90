!----------------------------------------------------------------------
!> @brief Flow field array and physical simualtion parameters
!----------------------------------------------------------------------
module flow
use velocityGrid, only: PI
implicit none
save

! flow parameters, to be read from NML: flowNml
double precision :: Kn, pressDrop, accom

double precision :: mu

! double precision, parameter :: Kn = 1.0d0
! double precision, parameter :: mu = dsqrt(PI/2.d0)/Kn
! double precision, parameter :: PressDrop=1.0d-1
! double precision, parameter :: accom = 1.d0

double precision, DIMENSION(:,:,:), ALLOCATABLE :: f
double precision, DIMENSION(:), ALLOCATABLE :: Rho, Ux, Uy, Uz
double precision :: mass

contains
    subroutine setupFlow     
        use physicalGrid
        use velocityGrid
        use mpiParams, only: f_west_snd,  f_east_snd,  f_west_rcv,  f_east_rcv, &
                             f_suth_snd,  f_noth_snd,  f_suth_rcv,  f_noth_rcv, &
                             f_back_snd,  f_frnt_snd,  f_back_rcv,  f_frnt_rcv
        implicit none
        integer :: i, j, k, l, ll

        ALLOCATE(f(Ntotal,Nc8,1:8))

        ALLOCATE(Rho(Ntotal), Ux(Ntotal), Uy(Ntotal), Uz(Ntotal))

        mu = dsqrt(PI/2.d0)/Kn

        f=0.d0
        Rho = 0.d0
        Ux = 0.d0
        Uy = 0.d0
        Uz = 0.d0
        mass = 1.d0


        Do k=zlg,zug
            Do j=ylg,yug
                Do i=xlg,xug
                    l=(i-xlg+1)+(j-ylg)*Nxtotal+(k-zlg)*Nxytotal
                    if (image(i,j,k)/=solid) then
            !           Rho(l)=pressDrop*(i/2.d0-Nx)/Nx
                        Rho(l)= 0.d0
                    end if
                    do ll=1,8
                        f(l,:,ll)=w(:)*Rho(l) !Check
                    enddo
                 Enddo
            Enddo
        Enddo
    end subroutine setupFlow

    subroutine allocateBuf
        use physicalGrid
        use velocityGrid
        use mpiParams
        implicit none

        ! Now we calculate the snd/rcv buffer sizes then we can allocate them
        ! Send buffer contains f* of first two real layers AND ALSO f*w of 
        ! 3-fold corners that inside the first two layers
        westSndSize = Nytotal*Nztotal*ghostLayers*Nc + westN3corner_snd*Nc
        westRcvSize = Nytotal*Nztotal*ghostLayers*Nc + westN3corner_rcv*Nc
        eastSndSize = Nytotal*Nztotal*ghostLayers*Nc + eastN3corner_snd*Nc
        eastRcvSize = Nytotal*Nztotal*ghostLayers*Nc + eastN3corner_rcv*Nc
        suthSndSize = Nxtotal*Nztotal*ghostLayers*Nc + suthN3corner_snd*Nc
        suthRcvSize = Nxtotal*Nztotal*ghostLayers*Nc + suthN3corner_rcv*Nc
        nothSndSize = Nxtotal*Nztotal*ghostLayers*Nc + nothN3corner_snd*Nc
        nothRcvSize = Nxtotal*Nztotal*ghostLayers*Nc + nothN3corner_rcv*Nc
        backSndSize = Nxtotal*Nytotal*ghostLayers*Nc + backN3corner_snd*Nc
        backRcvSize = Nxtotal*Nytotal*ghostLayers*Nc + backN3corner_rcv*Nc
        frntSndSize = Nxtotal*Nytotal*ghostLayers*Nc + frntN3corner_snd*Nc
        frntRcvSize = Nxtotal*Nytotal*ghostLayers*Nc + frntN3corner_rcv*Nc
        ! Allocate the buffer
        allocate( f_west_snd(westSndSize) )
        allocate( f_west_rcv(westRcvSize) )
        allocate( f_east_snd(eastSndSize) )
        allocate( f_east_rcv(eastRcvSize) )
        allocate( f_suth_snd(suthSndSize) )
        allocate( f_suth_rcv(suthRcvSize) )
        allocate( f_noth_snd(nothSndSize) )
        allocate( f_noth_rcv(nothRcvSize) )
        allocate( f_back_snd(backSndSize) )
        allocate( f_back_rcv(backRcvSize) )
        allocate( f_frnt_snd(frntSndSize) )
        allocate( f_frnt_rcv(frntRcvSize) )
        ! Allways initialize allocated arrays
        f_west_rcv = 0.d0
        f_east_rcv = 0.d0
        f_noth_rcv = 0.d0
        f_suth_rcv = 0.d0
        f_back_rcv = 0.d0
        f_frnt_rcv = 0.d0
        f_west_snd = 0.d0
        f_east_snd = 0.d0
        f_noth_snd = 0.d0
        f_suth_snd = 0.d0
        f_back_snd = 0.d0
        f_frnt_snd = 0.d0
    end subroutine allocateBuf
end module flow
