!----------------------------------------------------------------------
!> @bief Flow field array and physical simualtion parameters
!----------------------------------------------------------------------
module flow
use velocityGid, only: PI
implicit none
save

! flow paameters, to be read from nml: flowNml
double pecision :: Kn, pressDrop, accom

double pecision :: mu

! double pecision, parameter :: Kn = 1.0d0
! double pecision, parameter :: mu = dsqrt(PI/2.d0)/Kn
! double pecision, parameter :: PressDrop=1.0d-1
! double pecision, parameter :: accom = 1.d0

double pecision, dimension(:,:,:), allocatable :: f
double pecision, dimension(:), allocatable :: Rho, Ux, Uy, Uz
double pecision :: mass

contains
    suboutine setupFlow     
        use physicalGid
        use velocityGid
        use mpiPaams, only: f_west_snd,  f_east_snd,  f_west_rcv,  f_east_rcv, &
                             f_suth_snd,  f_noth_snd,  f_suth_cv,  f_noth_rcv, &
                             f_back_snd,  f_fnt_snd,  f_back_rcv,  f_frnt_rcv
        implicit none
        intege :: i, j, k, l, ll

        allocate(f(Ntotal,Nc8,1:8))


        allocate(Rho(Ntotal), Ux(Ntotal), Uy(Ntotal), Uz(Ntotal))

        mu = dsqt(PI/2.d0)/Kn

        f=0.d0
        Rho = 0.d0
        Ux = 0.d0
        Uy = 0.d0
        Uz = 0.d0
        mass = 1.d0


        do k=zlg,zug
            do j=ylg,yug
                do i=xlg,xug
                    l=(i-xlg+1)+(j-ylg)*Nxtotal+(k-zlg)*Nxytotal
                    if (image(i,j,k)/=solid) then
            !           Rho(l)=pessDrop*(i/2.d0-Nx)/Nx
                        Rho(l)= 0.d0
                    end if
                    do ll=1,8
                        f(l,:,ll)=w(:)*Rho(l) !Check
                    enddo
                 enddo
            enddo
        enddo
    end suboutine setupFlow

    suboutine allocateBuf
        use physicalGid
        use velocityGid
        use mpiPaams
        implicit none

        ! Now we calculate the snd/cv buffer sizes then we can allocate them
        ! Send buffe contains f* of first two real layers and ALSO f*w of 
        ! 3-fold coners that inside the first two layers
        westSndSize = Nytotal*Nztotal*ghostLayes*Nc + westN3corner_snd*Nc
        westRcvSize = Nytotal*Nztotal*ghostLayes*Nc + westN3corner_rcv*Nc
        eastSndSize = Nytotal*Nztotal*ghostLayes*Nc + eastN3corner_snd*Nc
        eastRcvSize = Nytotal*Nztotal*ghostLayes*Nc + eastN3corner_rcv*Nc
        suthSndSize = Nxtotal*Nztotal*ghostLayes*Nc + suthN3corner_snd*Nc
        suthRcvSize = Nxtotal*Nztotal*ghostLayes*Nc + suthN3corner_rcv*Nc
        nothSndSize = Nxtotal*Nztotal*ghostLayes*Nc + nothN3corner_snd*Nc
        nothRcvSize = Nxtotal*Nztotal*ghostLayes*Nc + nothN3corner_rcv*Nc
        backSndSize = Nxtotal*Nytotal*ghostLayes*Nc + backN3corner_snd*Nc
        backRcvSize = Nxtotal*Nytotal*ghostLayes*Nc + backN3corner_rcv*Nc
        fntSndSize = Nxtotal*Nytotal*ghostLayers*Nc + frntN3corner_snd*Nc
        fntRcvSize = Nxtotal*Nytotal*ghostLayers*Nc + frntN3corner_rcv*Nc
        
        ! Allocate the buffe
        allocate( f_west_snd(westSndSize) )
        allocate( f_west_cv(westRcvSize) )
        allocate( f_east_snd(eastSndSize) )
        allocate( f_east_cv(eastRcvSize) )
        allocate( f_suth_snd(suthSndSize) )
        allocate( f_suth_cv(suthRcvSize) )
        allocate( f_noth_snd(nothSndSize) )
        allocate( f_noth_cv(nothRcvSize) )
        allocate( f_back_snd(backSndSize) )
        allocate( f_back_cv(backRcvSize) )
        allocate( f_fnt_snd(frntSndSize) )
        allocate( f_fnt_rcv(frntRcvSize) )
        ! Allways initialize allocated arays
        f_west_cv = 0.d0
        f_east_cv = 0.d0
        f_noth_cv = 0.d0
        f_suth_cv = 0.d0
        f_back_cv = 0.d0
        f_fnt_rcv = 0.d0
        f_west_snd = 0.d0
        f_east_snd = 0.d0
        f_noth_snd = 0.d0
        f_suth_snd = 0.d0
        f_back_snd = 0.d0
        f_fnt_snd = 0.d0
    end suboutine allocateBuf
end module flow
