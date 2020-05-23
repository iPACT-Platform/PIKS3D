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
! module    : flow
!-------------------------------------------------------------------------------
! This is a module for the flow parameters of 3D DVM parallel solver. 
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Comput. Phys. Commun., 234 (2019), pp. 14-25
!
!	Initial conditions for velocity distribution function, velocity are set.
!-------------------------------------------------------------------------------

module flow
use velocityGrid, only: PI
implicit none
save

! flow parameters, to be read from NML: flowNml
! Kn : Knudsen number defined in Eq.(1) of [1]
! pressDrop : pressure drop applied at inlet/outlet Eq.(10) of [1]
! accom : tangential momentum accommodation coefficient (TMAC)  Eq.(8) of [1]
double precision :: Kn, pressDrop, accom
! mu : coefficient associated with collision frequency in RHS of Eq.(4) of [1]
double precision :: mu

! double precision, parameter :: Kn = 1.0d0
! double precision, parameter :: mu = dsqrt(PI/2.d0)/Kn
! double precision, parameter :: PressDrop=1.0d-1
! double precision, parameter :: accom = 1.d0

! f(spatial_id,velocity_id,sweeppath_id) : velocity distribution function in Eq.(5) of [1]
double precision, DIMENSION(:,:,:), ALLOCATABLE :: f
! Rho : number density in Eq.(6) of [1]
! Ux, Uy, Uz : three component of velocity vector U in Eq.(6) of [1]
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
