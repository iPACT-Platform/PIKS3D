module solver
    use flow
    use velocityGrid
    use physicalGrid
    use mpiParams

    implicit none

    ! to be read from NML: solverNml
    double precision :: eps
    integer :: maxStep
    integer :: chkConvergeStep
    integer :: saveStep
    logical :: saveLast
    integer :: saveFormat ! 1 (default) for vti, 2 for tecplot, 3 for vtk

    integer :: iStep
    double precision :: error
    double precision :: permeability

    contains
    subroutine iterate
        implicit none
        include "mpif.h"

        integer :: i, j, k, l, ii, jj, kk, m
        integer :: localidMsnd, localidPsnd, localidMrcv, localidPrcv, packid
        INTEGER :: MPI_ERR
        INTEGER :: MPI_REQ_X(4), MPI_REQ_Y(4), MPI_REQ_Z(4)
        INTEGER :: MPI_STAT(MPI_STATUS_SIZE,6)
        integer :: xfsize, yfsize, zfsize, locB
        double precision :: fEq, RhoWall, RhoWall2, RhoWall3
        double precision, dimension(Nc8,1:8) :: fwZ
        
        ! for mpi error handling
        integer :: errMsgLen, errTmp
        character(len=MPI_MAX_ERROR_STRING) :: errMsg 

        ! The size for storing two layers of nodes in buffer 
        xfsize = Nytotal*Nztotal*Nc*ghostLayers
        yfsize = Nxtotal*Nztotal*Nc*ghostLayers
        zfsize = Nxtotal*Nytotal*Nc*ghostLayers

        MPI_REQ_X = MPI_REQUEST_NULL
        MPI_REQ_Y = MPI_REQUEST_NULL
        MPI_REQ_Z = MPI_REQUEST_NULL

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(l, i, j, k,m, ii, jj, kk, fEq, RhoWall, RhoWall2, RhoWall3)&
!$OMP PRIVATE(fwZ)&
!$OMP PRIVATE(localidPsnd, localidPrcv, localidMsnd, localidMrcv, packid)


!$OMP SINGLE
        ! Start Recieving
        CALL MPI_IRECV( f_east_rcv, eastRcvSize, MPI_DOUBLE_PRECISION, east, TAG1, &
                        MPI_COMM_VGRID, MPI_REQ_X(1), MPI_ERR )
        CALL MPI_IRECV( f_west_rcv, westRcvSize, MPI_DOUBLE_PRECISION, west, TAG2, &
                        MPI_COMM_VGRID, MPI_REQ_X(2), MPI_ERR )
        CALL MPI_IRECV( f_noth_rcv, nothRcvSize, MPI_DOUBLE_PRECISION, noth, TAG3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(1), MPI_ERR )
        CALL MPI_IRECV( f_suth_rcv, suthRcvSize, MPI_DOUBLE_PRECISION, suth, TAG4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(2), MPI_ERR )
        CALL MPI_IRECV( f_frnt_rcv, frntRcvSize, MPI_DOUBLE_PRECISION, frnt, TAG5, &
                        MPI_COMM_VGRID, MPI_REQ_Z(1), MPI_ERR )
        CALL MPI_IRECV( f_back_rcv, backRcvSize, MPI_DOUBLE_PRECISION, back, TAG6, &
                        MPI_COMM_VGRID, MPI_REQ_Z(2), MPI_ERR )

        CALL MPI_ISEND( f_west_snd, westSndSize, MPI_DOUBLE_PRECISION, west, TAG1, &
                        MPI_COMM_VGRID, MPI_REQ_X(3), MPI_ERR )
        CALL MPI_ISEND( f_east_snd, eastSndSize, MPI_DOUBLE_PRECISION, east, TAG2, &
                        MPI_COMM_VGRID, MPI_REQ_X(4), MPI_ERR )
        CALL MPI_ISEND( f_suth_snd, suthSndSize, MPI_DOUBLE_PRECISION, suth, TAG3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(3), MPI_ERR )
        CALL MPI_ISEND( f_noth_snd, nothSndSize, MPI_DOUBLE_PRECISION, noth, TAG4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(4), MPI_ERR )
        CALL MPI_ISEND( f_back_snd, backSndSize, MPI_DOUBLE_PRECISION, back, TAG5, &
                        MPI_COMM_VGRID, MPI_REQ_Z(3), MPI_ERR )
        CALL MPI_ISEND( f_frnt_snd, frntSndSize, MPI_DOUBLE_PRECISION, frnt, TAG6, &
                        MPI_COMM_VGRID, MPI_REQ_Z(4), MPI_ERR )
!$OMP END SINGLE NOWAIT

!------------------------------------------------------------------------
!           In the 1st group of direction cx>0 & cy>0 & cz>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC) 
    Do l=1,Nc8
        Do i=1,Nstencil1
            !k=dir1(i)
            ii = dir1(1,i)
            jj = dir1(2,i)
            kk = dir1(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction x-y
            if (image(ii,jj-1,kk)==wallEN) then
                f(k-Nxtotal,l,1)=f(k-Nxtotal,l,3)
            elseif ((image(ii,jj-1,kk)==wallENF) &
                .or.(image(ii,jj-1,kk)==wallENB)) then
                f(k-Nxtotal,l,1)=fw(which_corner(ii,jj-1,kk),l,1)
            endif

            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk-1)==wallNF) then
                f(k-Nxytotal,l,1)=f(k-Nxytotal,l,8)
            elseif (image(ii,jj,kk-1)==wallEF) then
                f(k-Nxytotal,l,1)=f(k-Nxytotal,l,6)
            elseif ((image(ii,jj,kk-1)==wallWNF) &
               .or. (image(ii,jj,kk-1)==wallESF)) then
                f(k-Nxytotal,l,1)=fw(which_corner(ii,jj,kk-1),l,1)
            else if (image(ii,jj,kk-1)==wallENF) then
                f(k-Nxytotal,l,1)=f(k-Nxytotal,l,7)
            endif

            fEq=w(l)*(Rho(k)+1.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)+cz(l)*Uz(k)))
            f(k,l,1)=(mu*(fEq-0.5d0*f(k,l,1)) &
            &        + cx(l)*coef1(i,2)*f(k-1,l,1) &
            &        + cx(l)*coef1(i,3)*f(k-2,l,1) &
            &        + cy(l)*coef1(i,5)*f(k-Nxtotal,l,1) &
            &        + cy(l)*coef1(i,6)*f(k-2*Nxtotal,l,1) &
            &        + cz(l)*coef1(i,8)*f(k-Nxytotal,l,1) &
            &        + cz(l)*coef1(i,9)*f(k-2*Nxytotal,l,1) &
            & )/(0.5d0*mu+cx(l)*coef1(i,1)+cy(l)*coef1(i,4)+cz(l)*coef1(i,7))
        End do
    End do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 2nd group of direction cx<0 & cy>0 & cz>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    Do l=1,Nc8
        Do i=1,Nstencil2
            !k=dir1(i)
            ii = dir2(1,i)
            jj = dir2(2,i)
            kk = dir2(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected x&y-direction
            if (image(ii,jj-1,kk)==wallWN) then
                f(k-Nxtotal,l,2)=f(k-Nxtotal, l,4)
            elseif ((image(ii,jj-1,kk)==wallWNF) &
               .or. (image(ii,jj-1,kk)==wallWNB)) then
                f(k-Nxtotal,l,2)=fw(which_corner(ii,jj-1,kk),l,2)
            endif
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk-1)==wallNF) then
                f(k-Nxytotal,l,2)=f(k-Nxytotal,l,7)
            elseif (image(ii,jj,kk-1)==wallWF) then
                f(k-Nxytotal,l,2)=f(k-Nxytotal,l,5)
            elseif ((image(ii,jj,kk-1)==wallENF) &
               .or. (image(ii,jj,kk-1)==wallWSF)) then
                f(k-Nxytotal,l,2)=fw(which_corner(ii,jj,kk-1),l,2)
            elseif (image(ii,jj,kk-1)==wallWNF) then
                f(k-Nxytotal,l,2)=f(k-Nxytotal,l,8)
            end if

            fEq=w(l)*(Rho(k)+1.d0*(-cx(l)*Ux(k)+cy(l)*Uy(k)+cz(l)*Uz(k)))
            f(k,l,2)=(mu*(fEq-0.5d0*f(k,l,2)) &
            &        - cx(l)*coef2(i,2)*f(k+1,l,2) &
            &        - cx(l)*coef2(i,3)*f(k+2,l,2) &
            &        + cy(l)*coef2(i,5)*f(k-Nxtotal,l,2) &
            &        + cy(l)*coef2(i,6)*f(k-2*Nxtotal,l,2) &
            &        + cz(l)*coef2(i,8)*f(k-Nxytotal,l,2) &
            &        + cz(l)*coef2(i,9)*f(k-2*Nxytotal,l,2) &
            & )/(0.5d0*mu-cx(l)*coef2(i,1)+cy(l)*coef2(i,4)+cz(l)*coef2(i,7))
        End do
    End do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 3rd group of direction cx<0 & cy<0 & cz>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    Do l=1,Nc8
        Do i=1,Nstencil3
            ii = dir3(1,i)
            jj = dir3(2,i)
            kk = dir3(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction
            if (image(ii,jj+1,kk)==wallWS) then
                f(k+Nxtotal,l,3)=f(k+Nxtotal,l,1)
            elseif ((image(ii,jj+1,kk)==wallWSF) &
               .or. (image(ii,jj+1,kk)==wallWSB)) then
                f(k+Nxtotal,l,3)=fw(which_corner(ii,jj+1,kk),l,3)
            end if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk-1)==wallSF) then
                f(k-Nxytotal,l,3)=f(k-Nxytotal,l,6) 
            elseif (image(ii,jj,kk-1)==wallWF) then
                f(k-Nxytotal,l,3)=f(k-Nxytotal,l,8) 
            elseif ((image(ii,jj,kk-1)==wallWNF) &
               .or. (image(ii,jj,kk-1)==wallESF)) then
                f(k-Nxytotal,l,3)=fw(which_corner(ii,jj,kk-1),l,3)
            elseif (image(ii,jj,kk-1)==wallWSF) then
                f(k-Nxytotal,l,3)=f(k-Nxytotal,l,5)
            endif

            fEq=w(l)*(Rho(k)+1.d0*(-cx(l)*Ux(k)-cy(l)*Uy(k)+cz(l)*Uz(k)))
            f(k,l,3)=(mu*(fEq-0.5d0*f(k,l,3)) &
            &        - cx(l)*coef3(i,2)*f(k+1,l,3) &
            &        - cx(l)*coef3(i,3)*f(k+2,l,3) &
            &        - cy(l)*coef3(i,5)*f(k+Nxtotal,l,3) &
            &        - cy(l)*coef3(i,6)*f(k+2*Nxtotal,l,3) &
            &        + cz(l)*coef3(i,8)*f(k-Nxytotal,l,3) &
            &        + cz(l)*coef3(i,9)*f(k-2*Nxytotal,l,3) &
            & )/(0.5d0*mu-cx(l)*coef3(i,1)-cy(l)*coef3(i,4)+cz(l)*coef3(i,7))
        End do
    End do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 4th group of direction cx>0 & cy<0 & cz>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    Do l=1,Nc8
        Do i=1,Nstencil4
            ii = dir4(1,i)
            jj = dir4(2,i)
            kk = dir4(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction
            if (image(ii,jj+1,kk)==wallES) then
                f(k+Nxtotal,l,4)=f(k+Nxtotal,l,2)
            elseif ((image(ii,jj+1,kk)==wallESF) &
               .or. (image(ii,jj+1,kk)==wallESB)) then
                f(k+Nxtotal,l,4)=fw(which_corner(ii,jj+1,kk),l,4)
            endif
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk-1)==wallSF) then
                f(k-Nxytotal,l,4)=f(k-Nxytotal,l,5)
            elseif (image(ii,jj,kk-1)==wallEF) then
                f(k-Nxytotal,l,4)=f(k-Nxytotal,l,7)
            elseif ((image(ii,jj,kk-1)==wallENF) &
                .or.(image(ii,jj,kk-1)==wallWSF)) then
                f(k-Nxytotal,l,4)=fw(which_corner(ii,jj,kk-1),l,4)
            elseif (image(ii,jj,kk-1)==wallESF) then
                f(k-Nxytotal,l,4)=f(k-Nxytotal,l,6) 
            endif

            fEq=w(l)*(Rho(k)+1.d0*(cx(l)*Ux(k)-cy(l)*Uy(k)+cz(l)*Uz(k)))
            f(k,l,4)=(mu*(fEq-0.5d0*f(k,l,4)) &
            &        + cx(l)*coef4(i,2)*f(k-1,l,4) &
            &        + cx(l)*coef4(i,3)*f(k-2,l,4) &
            &        - cy(l)*coef4(i,5)*f(k+Nxtotal,l,4) &
            &        - cy(l)*coef4(i,6)*f(k+2*Nxtotal,l,4) &
            &        + cz(l)*coef4(i,8)*f(k-Nxytotal,l,4) &
            &        + cz(l)*coef4(i,9)*f(k-2*Nxytotal,l,4) &
            & )/(0.5d0*mu+cx(l)*coef4(i,1)-cy(l)*coef4(i,4)+cz(l)*coef4(i,7))
        End do
    End do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 5th group of direction cx>0 & cy>0 & cz<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    Do l=1,Nc8
        Do i=1,Nstencil5
            ii = dir5(1,i)
            jj = dir5(2,i)
            kk = dir5(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction x-y
            if (image(ii,jj-1,kk)==wallEN) then
                f(k-Nxtotal,l,5)=f(k-Nxtotal,l,7)
            elseif ((image(ii,jj-1,kk)==wallENF) &
               .or. (image(ii,jj-1,kk)==wallENB)) then
                f(k-Nxtotal,l,5)=fw(which_corner(ii,jj-1,kk),l,5)
            endif
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk+1)==wallNB) then
                f(k+Nxytotal,l,5)=f(k+Nxytotal,l,4)
            elseif (image(ii,jj,kk+1)==wallEB) then
                f(k+Nxytotal,l,5)=f(k+Nxytotal,l,2)
            elseif ((image(ii,jj,kk+1)==wallESB) &
               .or. (image(ii,jj,kk+1)==wallWNB)) then
                f(k+Nxytotal,l,5)=fw(which_corner(ii,jj,kk+1),l,5)
            elseif (image(ii,jj,kk+1)==wallENB) then
                f(k+Nxytotal,l,5)=f(k+Nxytotal,l,3)
            endif

            fEq=w(l)*(Rho(k)+1.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)-cz(l)*Uz(k)))
            f(k,l,5)=(mu*(fEq-0.5d0*f(k,l,5)) &
            &        + cx(l)*coef5(i,2)*f(k-1,l,5) &
            &        + cx(l)*coef5(i,3)*f(k-2,l,5) &
            &        + cy(l)*coef5(i,5)*f(k-Nxtotal,l,5) &
            &        + cy(l)*coef5(i,6)*f(k-2*Nxtotal,l,5) &
            &        - cz(l)*coef5(i,8)*f(k+Nxytotal,l,5) &
            &        - cz(l)*coef5(i,9)*f(k+2*Nxytotal,l,5) &
            & )/(0.5d0*mu+cx(l)*coef5(i,1)+cy(l)*coef5(i,4)-cz(l)*coef5(i,7))
        End do
    End do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 6th group of direction cx<0 & cy>0 & cz<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    Do l=1,Nc8
        Do i=1,Nstencil6
            ii = dir6(1,i)
            jj = dir6(2,i)
            kk = dir6(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected x&y-direction
            if (image(ii,jj-1,kk)==wallWN) then
                f(k-Nxtotal,l,6)=f(k-Nxtotal,l,8)
            elseif ((image(ii,jj-1,kk)==wallWNF) &
                .or.(image(ii,jj-1,kk)==wallWNB)) then
                f(k-Nxtotal,l,6)=fw(which_corner(ii,jj-1,kk),l,6)
            endif
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk+1)==wallNB) then
                f(k+Nxytotal,l,6)=f(k+Nxytotal,l,3)
            elseif (image(ii,jj,kk+1)==wallWB) then
                f(k+Nxytotal,l,6)=f(k+Nxytotal,l,1)
            elseif ((image(ii,jj,kk+1)==wallENB) &
               .or. (image(ii,jj,kk+1)==wallWSB)) then
                f(k+Nxytotal,l,6)=fw(which_corner(ii,jj,kk+1),l,6)
            elseif (image(ii,jj,kk+1)==wallWNB) then
                f(k+Nxytotal,l,6)=f(k+Nxytotal,l,4)
            endif

            fEq=w(l)*(Rho(k)+1.d0*(-cx(l)*Ux(k)+cy(l)*Uy(k)-cz(l)*Uz(k)))
            f(k,l,6)=(mu*(fEq-0.5d0*f(k,l,6)) &
            &        - cx(l)*coef6(i,2)*f(k+1,l,6) &
            &        - cx(l)*coef6(i,3)*f(k+2,l,6) &
            &        + cy(l)*coef6(i,5)*f(k-Nxtotal,l,6) &
            &        + cy(l)*coef6(i,6)*f(k-2*Nxtotal,l,6) &
            &        - cz(l)*coef6(i,8)*f(k+Nxytotal,l,6) &
            &        - cz(l)*coef6(i,9)*f(k+2*Nxytotal,l,6) &
            & )/(0.5d0*mu-cx(l)*coef6(i,1)+cy(l)*coef6(i,4)-cz(l)*coef6(i,7))
        End do
    End do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 7th group of direction cx<0 & cy<0 & cz<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    Do l=1,Nc8
        Do i=1,Nstencil7
            ii = dir7(1,i)
            jj = dir7(2,i)
            kk = dir7(3,i) 
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction
            if (image(ii,jj+1,kk)==wallWS) then
                f(k+Nxtotal,l,7)=f(k+Nxtotal,l,8)
            elseif ((image(ii,jj+1,kk)==wallWSF) &
               .or. (image(ii,jj+1,kk)==wallWSB)) then
                f(k+Nxtotal,l,7)=fw(which_corner(ii,jj+1,kk),l,7)
            endif
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk+1)==wallSB) then
                f(k+Nxytotal,l,7)=f(k+Nxytotal,l,2)
            elseif (image(ii,jj,kk+1)==wallWB) then
                f(k+Nxytotal,l,7)=f(k+Nxytotal,l,4)
            elseif ((image(ii,jj,kk+1)==wallWNB) &
               .or. (image(ii,jj,kk+1)==wallESB)) then
                f(k+Nxytotal,l,7)=fw(which_corner(ii,jj,kk+1),l,7)
            else if (image(ii,jj,kk+1)==wallWSB) then
                f(k+Nxytotal,l,7)=f(k+Nxytotal,l,1)
            endif

            fEq=w(l)*(Rho(k)+1.d0*(-cx(l)*Ux(k)-cy(l)*Uy(k)-cz(l)*Uz(k)))
            f(k,l,7)=(mu*(fEq-0.5d0*f(k,l,7)) &
            &        - cx(l)*coef7(i,2)*f(k+1,l,7) &
            &        - cx(l)*coef7(i,3)*f(k+2,l,7) &
            &        - cy(l)*coef7(i,5)*f(k+Nxtotal,l,7) &
            &        - cy(l)*coef7(i,6)*f(k+2*Nxtotal,l,7) &
            &        - cz(l)*coef7(i,8)*f(k+Nxytotal,l,7) &
            &        - cz(l)*coef7(i,9)*f(k+2*Nxytotal,l,7) &
            & )/(0.5d0*mu-cx(l)*coef7(i,1)-cy(l)*coef7(i,4)-cz(l)*coef7(i,7))
        End do
    End do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 8th group of direction cx>0 & cy<0 & cz<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    Do l=1,Nc8
        Do i=1,Nstencil8
            ii = dir8(1,i)
            jj = dir8(2,i)
            kk = dir8(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction
            if (image(ii,jj+1,kk)==wallES) then
                f(k+Nxtotal,l,8)=f(k+Nxtotal,l,6)
            elseif ((image(ii,jj+1,kk)==wallESF) &
               .or. (image(ii,jj+1,kk)==wallESB)) then
                f(k+Nxtotal,l,8)=fw(which_corner(ii,jj+1,kk),l,8)
            endif
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk+1)==wallSB) then
                f(k+Nxytotal,l,8)=f(k+Nxytotal,l,1)
            elseif (image(ii,jj,kk+1)==wallEB) then
                f(k+Nxytotal,l,8)=f(k+Nxytotal,l,3)
            elseif ((image(ii,jj,kk+1)==wallENB) &
               .or. (image(ii,jj,kk+1)==wallWSB)) then
                f(k+Nxytotal,l,8)=fw(which_corner(ii,jj,kk+1),l,8)
            elseif (image(ii,jj,kk+1)==wallESB) then
                f(k+Nxytotal,l,8)=f(k+Nxytotal,l,2)
            endif

            fEq=w(l)*(Rho(k)+1.d0*(cx(l)*Ux(k)-cy(l)*Uy(k)-cz(l)*Uz(k)))
            f(k,l,8)=(mu*(fEq-0.5d0*f(k,l,8)) &
            &        + cx(l)*coef8(i,2)*f(k-1,l,8) &
            &        + cx(l)*coef8(i,3)*f(k-2,l,8) &
            &        - cy(l)*coef8(i,5)*f(k+Nxtotal,l,8) &
            &        - cy(l)*coef8(i,6)*f(k+2*Nxtotal,l,8) &
            &        - cz(l)*coef8(i,8)*f(k+Nxytotal,l,8) &
            &        - cz(l)*coef8(i,9)*f(k+2*Nxytotal,l,8) &
            & )/(0.5d0*mu+cx(l)*coef8(i,1)-cy(l)*coef8(i,4)-cz(l)*coef8(i,7))
        End do
    End do
!$OMP END DO 

!$OMP SINGLE
        ! Wait until send and recv done
        CALL MPI_WAITALL(4, MPI_REQ_X, MPI_STAT, MPI_ERR)
        CALL MPI_WAITALL(4, MPI_REQ_Y, MPI_STAT, MPI_ERR)
        CALL MPI_WAITALL(4, MPI_REQ_Z, MPI_STAT, MPI_ERR)
!$OMP END SINGLE 

!$OMP DO COLLAPSE(3)
    ! pack/unpack X dir buffer
    do k = 1,Nztotal
        do j = 1, Nytotal
            do i = 1, ghostLayers
                ! pack order, 
                !for dir: (2,3,6,7)
                localidMsnd = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i+ghostLayers
                localidMrcv = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i
                !for dir: (1,4,5,8)
                localidPsnd = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i+Nxsub
                localidPrcv = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i+ghostLayers+Nxsub
        
                packid = (k-1)*Nytotal*ghostLayers*Nc &
                       + (j-1)*ghostLayers*Nc &
                       + (i-1)*Nc
                do l = 1,8
                    f_west_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidMsnd,:,l)
                    f(localidMrcv,:,l) = f_west_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                    f_east_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidPsnd,:,l)
                    f(localidPrcv,:,l) = f_east_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                enddo !l
            enddo !i
        enddo !j
    enddo !k
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3)
    ! pack/unpack Y dir buffer
    do k = 1, Nztotal
        do j = 1, ghostLayers
            do i = 1, Nxtotal
                !for dir: (3,4,7,8)
                localidMsnd = (k-1)*Nxtotal*Nytotal + (j+ghostLayers-1)*Nxtotal + i
                localidMrcv = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i
                !for dir: (2,1,6,5)
                localidPsnd = (k-1)*Nxtotal*Nytotal + (j+Nysub-1)*Nxtotal + i
                localidPrcv = (k-1)*Nxtotal*Nytotal + (j+Nysub+ghostLayers-1)*Nxtotal + i
                packid = (k-1)*ghostLayers*Nxtotal*Nc &
                       + (j-1)*Nxtotal*Nc &
                       + (i-1)*Nc
                do l = 1,8
                    f_suth_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidMsnd,:,l)
                    f(localidMrcv,:,l) = f_suth_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                    f_noth_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidPsnd,:,l)
                    f(localidPrcv,:,l) = f_noth_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                enddo !l
            enddo !i
        enddo !j
    enddo !k
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3)
    ! pack/unpack Z dir buffer
    do k = 1, ghostLayers
        do j = 1, Nytotal
            do i = 1, Nxtotal
                ! pack order, 
                !for dir: (5,6,7,8)
                localidMsnd = (k+ghostLayers-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                localidMrcv = (k-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                !for dir: (1,2,3,4)
                localidPsnd = (k+Nzsub-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                localidPrcv = (k+Nzsub+ghostLayers-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                packid = (k-1)*Nytotal*Nxtotal*Nc &
                       + (j-1)*Nxtotal*Nc &
                       + (i-1)*Nc
                do l = 1,8
                    f_back_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidMsnd,:,l)
                    f(localidMrcv,:,l) = f_back_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                    f_frnt_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidPsnd,:,l)
                    f(localidPrcv,:,l) = f_frnt_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                enddo !l                
            enddo !i
        enddo !j
    enddo !k
!$OMP END DO NOWAIT


!------------------------------------------------------------------------
! pack/unpack 3-fold corners at west/east
!------------------------------------------------------------------------
!$OMP SECTIONS
!$OMP SECTION
!!$OMP DO
    ! pack 3-fold corners send buffer at west
    do m = 1, westN3corner_snd
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            f_west_snd(locB+1:locB+Nc8) = fw(map3CorWsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! pack 3-fold corners send buffer at east
    do m = 1, eastN3corner_snd
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            f_east_snd(locB+1:locB+Nc8) = fw(map3CorEsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold corners rcv buffer at west
    do m = 1, westN3corner_rcv
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CorWrcv(m),:,l) = f_west_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold corners rcv buffer at east
    do m = 1, eastN3corner_rcv
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CorErcv(m),:,l) = f_east_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO NOWAIT

!------------------------------------------------------------------------
! pack/unpack 3-fold corners at noth/suth
!------------------------------------------------------------------------
!$OMP SECTION
!!$OMP DO
    ! pack 3-fold corners send buffer at suth
    do m = 1, suthN3corner_snd
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            f_suth_snd(locB+1:locB+Nc8) = fw(map3CorSsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO

!$OMP SECTION
!!$OMP DO
    ! pack 3-fold corners send buffer at noth
    do m = 1, nothN3corner_snd
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            f_noth_snd(locB+1:locB+Nc8) = fw(map3CorNsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold corners rcv buffer at suth
    do m = 1, suthN3corner_rcv
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CorSrcv(m),:,l) = f_suth_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold corners rcv buffer at noth
    do m = 1, nothN3corner_rcv
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CorNrcv(m),:,l) = f_noth_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO

!------------------------------------------------------------------------
! pack/unpack 3-fold corners at back/frnt
!------------------------------------------------------------------------
!$OMP SECTION
!!$OMP DO
    ! pack 3-fold corners send buffer at back
    do m = 1, backN3corner_snd
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            f_back_snd(locB+1:locB+Nc8) = fw(map3CorBsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! pack 3-fold corners send buffer at frnt
    do m = 1, frntN3corner_snd
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            f_frnt_snd(locB+1:locB+Nc8) = fw(map3CorFsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold corners rcv buffer at back
    do m = 1, backN3corner_rcv
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CorBrcv(m),:,l) = f_back_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold corners rcv buffer at frnt
    do m = 1, frntN3corner_rcv
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CorFrcv(m),:,l) = f_frnt_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO NOWAIT
!$OMP END SECTIONS

!=======================================================================
!     Boundary condition on the flat wall
!=======================================================================
!$OMP DO SCHEDULE(STATIC)
    Do i=1,nWall
        k=vecWall(i)
        kk = k/Nxytotal + zlg ! to be checked
        jj = (k-(kk-zlg)*Nxytotal)/Nxtotal + ylg
        ii = k-(kk-zlg)*Nxytotal-(jj-ylg)*Nxtotal + xlg -1      
        j=which_corner(ii,jj,kk)
        RhoWall=0.d0
        RhoWall2=0.d0
        RhoWall3=0.d0
        SELECT CASE (image(ii,jj,kk))
            CASE (wallE)
                Do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)
                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                Enddo
            CASE (wallW)
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)
                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                Enddo
            CASE (wallN)
                Do l=1,Nc8
                    f(k,l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    f(k,l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    f(k,l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    f(k,l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)
                    RhoWall=RhoWall+cy(l)*(f(k,l,3)+f(k,l,4)+f(k,l,7)+f(k,l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                Enddo
            CASE (wallS)
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    f(k,l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    f(k,l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    f(k,l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)
                    RhoWall=RhoWall+cy(l)*(f(k,l,1)+f(k,l,2)+f(k,l,5)+f(k,l,6))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                Enddo
            CASE (wallF)
                Do l=1,Nc8
                    f(k,l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    f(k,l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    f(k,l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)
                    RhoWall=RhoWall+cz(l)*(f(k,l,5)+f(k,l,6)+f(k,l,7)+f(k,l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                Enddo
            CASE (wallB)
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    f(k,l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    f(k,l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)
                    RhoWall=RhoWall+cz(l)*(f(k,l,1)+f(k,l,2)+f(k,l,3)+f(k,l,4))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                Enddo
!=======================================================================
!     Boundary condition on the 2-direction corner wall
!=======================================================================
!------------------------------------------------------------------------
!           No z-direction type wall
!------------------------------------------------------------------------
            CASE (wallEN)
                Do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)

                    fwZ(l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)

                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                    RhoWall2=RhoWall2+cy(l)*(fwZ(l,3)+fwZ(l,4)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write y
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    ! Store y to unused location
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                Enddo

            CASE (wallWN)
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)

                    fwZ(l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)

                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                    RhoWall2=RhoWall2+cy(l)*(fwZ(l,3)+fwZ(l,4)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8

                    ! Write x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Write y
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Store y to unused location
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                Enddo

            CASE (wallES)
                Do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    fwZ(l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)

                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                    RhoWall2=RhoWall2+cy(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,5)+fwZ(l,6))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write y
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    ! Store y to unused location
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                Enddo

            CASE (wallWS)
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    fwZ(l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)

                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                    RhoWall2=RhoWall2+cy(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,5)+fwZ(l,6))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Write y
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    ! Store y to unused location
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                Enddo

!------------------------------------------------------------------------
!           No y-direction type wall
!------------------------------------------------------------------------
            CASE (wallEF)
                Do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)

                    fwZ(l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)

                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                    RhoWall2=RhoWall2+cz(l)*(fwZ(l,5)+fwZ(l,6)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write z
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    ! Store z to unsed location
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                Enddo

            CASE (wallWF)
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)

                    fwZ(l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)

                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                    RhoWall2=RhoWall2+cz(l)*(fwZ(l,5)+fwZ(l,6)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Write z
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Store z to unsed location
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                Enddo

            CASE (wallEB)
                Do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    fwZ(l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)

                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                    RhoWall2=RhoWall2+cz(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,3)+fwZ(l,4))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write z
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    ! Store z to unsed location
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                Enddo

            CASE (wallWB)
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    fwZ(l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)

                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                    RhoWall2=RhoWall2+cz(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,3)+fwZ(l,4))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Write z
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    ! Store z
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                Enddo
!------------------------------------------------------------------------
!           No x-direction type wall
!------------------------------------------------------------------------
            CASE (wallNF)
                Do l=1,Nc8
                    f(k,l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    f(k,l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    f(k,l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    f(k,l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)

                    fwZ(l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)

                    RhoWall=RhoWall+cy(l)*(f(k,l,3)+f(k,l,4)+f(k,l,7)+f(k,l,8))
                    RhoWall2=RhoWall2+cz(l)*(fwZ(l,5)+fwZ(l,6)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write y
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write z
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Store z to unsed location
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                Enddo

            CASE (wallSF)
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    f(k,l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    f(k,l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    f(k,l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)

                    fwZ(l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)

                    RhoWall=RhoWall+cy(l)*(f(k,l,1)+f(k,l,2)+f(k,l,5)+f(k,l,6))
                    RhoWall2=RhoWall2+cz(l)*(fwZ(l,5)+fwZ(l,6)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write y
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    ! Write z
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    ! Store z to unsed location
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                Enddo

            CASE (wallNB)
                Do l=1,Nc8
                    f(k,l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    f(k,l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    f(k,l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    f(k,l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    fwZ(l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)

                    RhoWall=RhoWall+cy(l)*(f(k,l,3)+f(k,l,4)+f(k,l,7)+f(k,l,8))
                    RhoWall2=RhoWall2+cz(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,3)+fwZ(l,4))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write y
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write z
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    ! Store z to unsed location
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                Enddo

            CASE (wallSB)
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    f(k,l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    f(k,l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    f(k,l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    fwZ(l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)

                    RhoWall=RhoWall+cy(l)*(f(k,l,1)+f(k,l,2)+f(k,l,5)+f(k,l,6))
                    RhoWall2=RhoWall2+cz(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,3)+fwZ(l,4))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write y
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    ! Write z
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    ! Store z to unsed location
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                Enddo
!=======================================================================
!     Boundary condition on the 3-direction corner wall
!=======================================================================
            CASE (wallENF) !direction1
                Do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)

                    fw(j,l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    fw(j,l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    fw(j,l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    fw(j,l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)

                    fwZ(l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)

                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                    RhoWall2=RhoWall2+cy(l)*(fw(j,l,3)+fw(j,l,4)+fw(j,l,7)+fw(j,l,8))
                    RhoWall3=RhoWall3+cz(l)*(fwZ(l,5)+fwZ(l,6)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write y
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,3)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,7)
                    ! Write z
                    f(k,l,3)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,7)

                    ! Store y
                    fw(j,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,4)
                    fw(j,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,8)
                    ! Store z
                    fw(j,l,2)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    fw(j,l,4)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f(k,l,7)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,5)
                Enddo
            CASE (wallWNF) !direction2
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)

                    fw(j,l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    fw(j,l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    fw(j,l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    fw(j,l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)

                    fwZ(l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)

                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                    RhoWall2=RhoWall2+cy(l)*(fw(j,l,3)+fw(j,l,4)+fw(j,l,7)+fw(j,l,8))
                    RhoWall3=RhoWall3+cz(l)*(fwZ(l,5)+fwZ(l,6)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Write y
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,4)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,8)
                    ! Write z
                    f(k,l,4)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,8)

                    ! Store y
                    fw(j,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,3)
                    fw(j,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,7)
                    ! Store z
                    fw(j,l,1)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    fw(j,l,3)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f(k,l,8)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,6)
                Enddo
            CASE (wallWSF) !direction3
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)

                    fw(j,l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    fw(j,l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    fw(j,l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    fw(j,l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)

                    fwZ(l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)

                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                    RhoWall2=RhoWall2+cy(l)*(fw(j,l,1)+fw(j,l,2)+fw(j,l,5)+fw(j,l,6))
                    RhoWall3=RhoWall3+cz(l)*(fwZ(l,5)+fwZ(l,6)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Write y
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,1)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,5)
                    ! Write z
                    f(k,l,1)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,5)

                    ! Store y
                    fw(j,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,2)
                    fw(j,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,6)
                    ! Store z
                    fw(j,l,2)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    fw(j,l,4)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f(k,l,5)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,7)
                Enddo
            CASE (wallESF) !direction4
                Do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)

                    fw(j,l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    fw(j,l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    fw(j,l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    fw(j,l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)

                    fwZ(l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    fwZ(l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    fwZ(l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    fwZ(l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)

                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                    RhoWall2=RhoWall2+cy(l)*(fw(j,l,1)+fw(j,l,2)+fw(j,l,5)+fw(j,l,6))
                    RhoWall3=RhoWall3+cz(l)*(fwZ(l,5)+fwZ(l,6)+fwZ(l,7)+fwZ(l,8))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write y
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,2)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,6)
                    ! Write z
                    f(k,l,2)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,6)

                    ! Store y
                    fw(j,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,1)
                    fw(j,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,5)
                    ! Store z
                    fw(j,l,1)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    fw(j,l,3)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f(k,l,6)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,8)
                Enddo
            CASE (wallENB) !direction5
                Do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)

                    fw(j,l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    fw(j,l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    fw(j,l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    fw(j,l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    fwZ(l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)

                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                    RhoWall2=RhoWall2+cy(l)*(fw(j,l,3)+fw(j,l,4)+fw(j,l,7)+fw(j,l,8))
                    RhoWall3=RhoWall3+cz(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,3)+fwZ(l,4))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write y
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,3)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,7)
                    ! Write z
                    f(k,l,7)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,3)

                    ! Store y
                    fw(j,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,4)
                    fw(j,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,8)
                    ! Store z
                    fw(j,l,6)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    fw(j,l,8)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f(k,l,3)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,1)
                Enddo
            CASE (wallWNB) !direction6
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)

                    fw(j,l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    fw(j,l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    fw(j,l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    fw(j,l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    fwZ(l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)

                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                    RhoWall2=RhoWall2+cy(l)*(fw(j,l,3)+fw(j,l,4)+fw(j,l,7)+fw(j,l,8))
                    RhoWall3=RhoWall3+cz(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,3)+fwZ(l,4))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Write y
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,4)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,8)
                    ! Write z
                    f(k,l,8)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,4)

                    ! Store y
                    fw(j,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,3)
                    fw(j,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,7)
                    ! Store z
                    fw(j,l,5)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    fw(j,l,7)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f(k,l,4)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,2)
                Enddo
            CASE (wallWSB) !direction7
                Do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)

                    fw(j,l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    fw(j,l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    fw(j,l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    fw(j,l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    fwZ(l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)

                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                    RhoWall2=RhoWall2+cy(l)*(fw(j,l,1)+fw(j,l,2)+fw(j,l,5)+fw(j,l,6))
                    RhoWall3=RhoWall3+cz(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,3)+fwZ(l,4))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Write y
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,1)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,5)
                    ! Write z
                    f(k,l,5)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,1)

                    ! Store y
                    fw(j,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,2)
                    fw(j,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,6)
                    ! Store z
                    fw(j,l,6)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    fw(j,l,8)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f(k,l,1)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,3)
                Enddo
            CASE (wallESB) !direction8
                Do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)

                    fw(j,l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    fw(j,l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    fw(j,l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    fw(j,l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)

                    fwZ(l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    fwZ(l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    fwZ(l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    fwZ(l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)

                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                    RhoWall2=RhoWall2+cy(l)*(fw(j,l,1)+fw(j,l,2)+fw(j,l,5)+fw(j,l,6))
                    RhoWall3=RhoWall3+cz(l)*(fwZ(l,1)+fwZ(l,2)+fwZ(l,3)+fwZ(l,4))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Write y
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,2)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,6)
                    ! Write z
                    f(k,l,6)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,2)

                    ! Store y
                    fw(j,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,1)
                    fw(j,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,5)
                    ! Store z
                    fw(j,l,5)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    fw(j,l,7)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f(k,l,2)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,4)
                Enddo
        END SELECT
    End do
!$OMP END DO   


!--------------------------------------------------------
!> inlet/outlet
!--------------------------------------------------------
    if(xl==xmin) then ! inlet block (west most processor)
!$OMP DO
        Do k=zlg,zug
            Do j=ylg,yug
                i = xl
                l = (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                !inlet
                f(l,:,1)=f(l-1,:,1) + w(:)*PressDrop
                f(l,:,4)=f(l-1,:,4) + w(:)*PressDrop
                f(l,:,5)=f(l-1,:,5) + w(:)*PressDrop
                f(l,:,8)=f(l-1,:,8) + w(:)*PressDrop
            End do
        End do
!$OMP END DO
    endif

    if(xu==xmax) then ! outlet block (east most processor)
!$OMP DO
        Do k=zlg,zug
            Do j=ylg,yug
                i = xu
                l = (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                !outlet
                f(l,:,2)=f(l+1,:,2) - w(:)*PressDrop
                f(l,:,3)=f(l+1,:,3) - w(:)*PressDrop
                f(l,:,6)=f(l+1,:,6) - w(:)*PressDrop
                f(l,:,7)=f(l+1,:,7) - w(:)*PressDrop
            Enddo
        End do
!$OMP END DO            
    endif

!----------------------------------------------------
!> Symmetric BC
!----------------------------------------------------
    if(yl==ymin) then ! south sym BC
!$OMP DO
        Do k=zl,zu
            Do i=xl,xu
                l = (k-zlg)*Nxytotal + ghostLayers*Nxtotal + i-xlg+1
                f(l,:,2)=f(l,:,3)
                f(l,:,1)=f(l,:,4)
                f(l,:,6)=f(l,:,7)
                f(l,:,5)=f(l,:,8)
            enddo
        End do
!$OMP END DO 
    endif
    if(yu==ymax) then ! north sym BC
!$OMP DO
        Do k=zl,zu
            Do i=xl,xu
                l = (k-zlg)*Nxytotal + (ghostLayers+Nysub-1)*Nxtotal + i-xlg+1
                f(l,:,3)=f(l,:,2)
                f(l,:,4)=f(l,:,1)
                f(l,:,7)=f(l,:,6)
                f(l,:,8)=f(l,:,5)
            End do
        End do
!$OMP END DO
    endif
    if(zl==zmin) then ! back sym BC
!$OMP DO
        Do j=yl,yu
            Do i=xl,xu
                l = ghostLayers*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                f(l,:,1)=f(l,:,5)
                f(l,:,2)=f(l,:,6)
                f(l,:,3)=f(l,:,7)
                f(l,:,4)=f(l,:,8)
            enddo
        End do
!$OMP END DO 
    endif
    if(zu==zmax) then ! front sym BC
!$OMP DO
        Do j=yl,yu
            Do i=xl,xu
                l = (ghostLayers+Nzsub-1)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                f(l,:,5)=f(l,:,1)
                f(l,:,6)=f(l,:,2)
                f(l,:,7)=f(l,:,3)
                f(l,:,8)=f(l,:,4)
            End do
        End do
!$OMP END DO
    endif

!----------------------------------------------------
!> Update Macro
!----------------------------------------------------
!$OMP DO SCHEDULE(STATIC) 
    !Do k=1,Ntotal
    Do i=1,Nfluid
        k=mapF(i)
        Rho(k)=0.d0
        Ux(k)=0.d0
        Uy(k)=0.d0
        Uz(k)=0.d0
        Do l=1,Nc8
            Rho(k)=Rho(k)+f(k,l,1)+f(k,l,2)+f(k,l,3)+f(k,l,4)+f(k,l,5)+f(k,l,6)+f(k,l,7)+f(k,l,8)
            Ux(k)=Ux(k)+cx(l)*(f(k,l,1)-f(k,l,2)-f(k,l,3)+f(k,l,4)+f(k,l,5)-f(k,l,6)-f(k,l,7)+f(k,l,8))
            Uy(k)=Uy(k)+cy(l)*(f(k,l,1)+f(k,l,2)-f(k,l,3)-f(k,l,4)+f(k,l,5)+f(k,l,6)-f(k,l,7)-f(k,l,8))
            Uz(k)=Uz(k)+cz(l)*(f(k,l,1)+f(k,l,2)+f(k,l,3)+f(k,l,4)-f(k,l,5)-f(k,l,6)-f(k,l,7)-f(k,l,8))
        End do
    End do
!$OMP END DO 
!$OMP END PARALLEL  
    end subroutine iterate
end module solver
