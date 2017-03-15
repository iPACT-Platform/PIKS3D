!----------------------------------------------------------------------
!> @brief Flow field array and physical simualtion parameters
!----------------------------------------------------------------------
module flow
use velocityGrid, only: PI
implicit none
save

double precision, parameter :: Kn = 9d-5
double precision, parameter :: mu = dsqrt(PI)/2.0d0/Kn
double precision, parameter :: PressDrop=1.0d-3
double precision, parameter :: accom = 1.d0

double precision, DIMENSION(:,:), ALLOCATABLE :: f1,f2,f3,f4,f5,f6,f7,f8
double precision, DIMENSION(:), ALLOCATABLE :: Rho, Ux, Uy, Uz
double precision :: mass

contains
    subroutine setupFlow     
        use physicalGrid, only: Ntotal, Nxtotal, Nytotal, Nztotal, ghostLayers
        use velocityGrid, only: Nc, Nc8
        use mpiParams, only: f_west_snd,  f_east_snd,  f_west_rcv,  f_east_rcv, &
                             f_suth_snd,  f_noth_snd,  f_soth_rcv,  f_noth_rcv, &
                             f_back_snd,  f_frnt_snd,  f_back_rcv,  f_frnt_rdv
        implicit none

        ALLOCATE(f1(Ntotal,Nc8), f2(Ntotal,Nc8), f3(Ntotal,Nc8), f4(Ntotal,Nc8), &
                 f5(Ntotal,Nc8), f6(Ntotal,Nc8), f7(Ntotal,Nc8), f8(Ntotal,Nc8) )

        ALLOCATE(Rho(Ntotal), Ux(Ntotal), Uy(Ntotal), Uz(Ntotal))

        f1=0.d0
        Rho = 0.d0
        Ux = 0.d0
        Uy = 0.d0
        Uz = 0.d0

        mass = 1.d0

        ! Buffer for exchange data only half domain of the velocity
        allocate( f_west_snd(Nc/2*Nytotal*Nztotal*ghostLayers) ) 
        allocate( f_west_rcv(Nc/2*Nytotal*Nztotal*ghostLayers) ) 
        allocate( f_east_snd(Nc/2*Nytotal*Nztotal*ghostLayers) ) 
        allocate( f_east_rcv(Nc/2*Nytotal*Nztotal*ghostLayers) ) 

        allocate( f_noth_snd(Nc/2*Nxtotal*Nztotal*ghostLayers) ) 
        allocate( f_noth_rcv(Nc/2*Nxtotal*Nztotal*ghostLayers) )
        allocate( f_suth_snd(Nc/2*Nxtotal*Nztotal*ghostLayers) ) 
        allocate( f_suth_rcv(Nc/2*Nxtotal*Nztotal*ghostLayers) )

        allocate( f_back_snd(Nc/2*Nxtotal*Nytotal*ghostLayers) ) 
        allocate( f_back_rcv(Nc/2*Nxtotal*Nytotal*ghostLayers) ) 
        allocate( f_frnt_snd(Nc/2*Nxtotal*Nytotal*ghostLayers) ) 
        allocate( f_frnt_rcv(Nc/2*Nxtotal*Nytotal*ghostLayers) ) 

        ! Allways initialize allocated arrays
        f1_west_rcv = 0.d0
        f1_east_rcv = 0.d0
        f1_noth_rcv = 0.d0
        f1_suth_rcv = 0.d0
        f1_back_rcv = 0.d0
        f1_frnt_rcv = 0.d0
        f1_west_snd = 0.d0
        f1_east_snd = 0.d0
        f1_noth_snd = 0.d0
        f1_suth_snd = 0.d0
        f1_back_snd = 0.d0
        f1_frnt_snd = 0.d0

    end subroutine setupFlow
end module flow
