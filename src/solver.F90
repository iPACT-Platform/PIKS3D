module solve
    use flow
    use velocityGid
    use physicalGid
    use mpiPaams

    implicit none

    ! to be ead from nml: solverNml
    double pecision :: eps
    intege :: maxStep
    intege :: chkConvergeStep
    intege :: saveStep
    logical :: saveLast
    intege :: saveFormat ! 1 (default) for vti, 2 for tecplot, 3 for vtk

    intege :: iStep
    double pecision :: error
    double pecision :: permeability

    contains
    suboutine iterate
        implicit none
        include "mpif.h"

        intege :: i, j, k, l, ii, jj, kk, m
        intege :: localidMsnd, localidPsnd, localidMrcv, localidPrcv, packid
        intege :: MPI_ERR
        intege :: MPI_REQ_X(4), MPI_REQ_Y(4), MPI_REQ_Z(4)
        intege :: MPI_STAT(MPI_STATUS_SIZE,6)
        intege :: xfsize, yfsize, zfsize, locB
        double pecision :: fEq, RhoWall, RhoWall2, RhoWall3
        double pecision, dimension(Nc8,1:8) :: fwZ
        
        ! fo mpi error handling
        intege :: errMsgLen, errTmp
        chaacter(len=MPI_MAX_ERROR_STRING) :: errMsg 

        ! The size fo storing two layers of nodes in buffer 
        xfsize = Nytotal*Nztotal*Nc*ghostLayes
        yfsize = Nxtotal*Nztotal*Nc*ghostLayes
        zfsize = Nxtotal*Nytotal*Nc*ghostLayes

        MPI_REQ_X = MPI_REQUEST_NULL
        MPI_REQ_Y = MPI_REQUEST_NULL
        MPI_REQ_Z = MPI_REQUEST_NULL

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(l, i, j, k,m, ii, jj, kk, fEq, RhoWall, RhoWall2, RhoWall3)&
!$OMP PRIVATE(fwZ)&
!$OMP PRIVATE(localidPsnd, localidPcv, localidMsnd, localidMrcv, packid)


!$OMP SINGLE
        ! Stat Recieving
        call MPI_IRECV( f_east_cv, eastRcvSize, MPI_DOUBLE_PRECISION, east, TAG1, &
                        MPI_COMM_VGRID, MPI_REQ_X(1), MPI_ERR )
        call MPI_IRECV( f_west_cv, westRcvSize, MPI_DOUBLE_PRECISION, west, TAG2, &
                        MPI_COMM_VGRID, MPI_REQ_X(2), MPI_ERR )
        call MPI_IRECV( f_noth_cv, nothRcvSize, MPI_DOUBLE_PRECISION, noth, TAG3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(1), MPI_ERR )
        call MPI_IRECV( f_suth_cv, suthRcvSize, MPI_DOUBLE_PRECISION, suth, TAG4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(2), MPI_ERR )
        call MPI_IRECV( f_fnt_rcv, frntRcvSize, MPI_DOUBLE_PRECISION, frnt, TAG5, &
                        MPI_COMM_VGRID, MPI_REQ_Z(1), MPI_ERR )
        call MPI_IRECV( f_back_cv, backRcvSize, MPI_DOUBLE_PRECISION, back, TAG6, &
                        MPI_COMM_VGRID, MPI_REQ_Z(2), MPI_ERR )

        call MPI_ISEND( f_west_snd, westSndSize, MPI_DOUBLE_PRECISION, west, TAG1, &
                        MPI_COMM_VGRID, MPI_REQ_X(3), MPI_ERR )
        call MPI_ISEND( f_east_snd, eastSndSize, MPI_DOUBLE_PRECISION, east, TAG2, &
                        MPI_COMM_VGRID, MPI_REQ_X(4), MPI_ERR )
        call MPI_ISEND( f_suth_snd, suthSndSize, MPI_DOUBLE_PRECISION, suth, TAG3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(3), MPI_ERR )
        call MPI_ISEND( f_noth_snd, nothSndSize, MPI_DOUBLE_PRECISION, noth, TAG4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(4), MPI_ERR )
        call MPI_ISEND( f_back_snd, backSndSize, MPI_DOUBLE_PRECISION, back, TAG5, &
                        MPI_COMM_VGRID, MPI_REQ_Z(3), MPI_ERR )
        call MPI_ISEND( f_fnt_snd, frntSndSize, MPI_DOUBLE_PRECISION, frnt, TAG6, &
                        MPI_COMM_VGRID, MPI_REQ_Z(4), MPI_ERR )
!$OMP END SINGLE NOWAIT

!------------------------------------------------------------------------
!           In the 1st goup of direction cx>0 & cy>0 & cz>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC) 
    do l=1,Nc8
        do i=1,Nstencil1
            !k=di1(i)
            ii = di1(1,i)
            jj = di1(2,i)
            kk = di1(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch fom reflected in x-direction (default) to in y-direction
            ! only fo group of velocity overlapped by two reflected direction x-y
            if (image(ii,jj-1,kk)==wallEN) then
                f(k-Nxtotal,l,1)=f(k-Nxtotal,l,3)
            elseif ((image(ii,jj-1,kk)==wallENF) &
                .o.(image(ii,jj-1,kk)==wallENB)) then
                f(k-Nxtotal,l,1)=fw(which_coner(ii,jj-1,kk),l,1)
            endif

            ! Switch fom reflected in x/y-direction to in z-direction
            ! only fo group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk-1)==wallNF) then
                f(k-Nxytotal,l,1)=f(k-Nxytotal,l,8)
            elseif (image(ii,jj,kk-1)==wallEF) then
                f(k-Nxytotal,l,1)=f(k-Nxytotal,l,6)
            elseif ((image(ii,jj,kk-1)==wallWNF) &
               .o. (image(ii,jj,kk-1)==wallESF)) then
                f(k-Nxytotal,l,1)=fw(which_coner(ii,jj,kk-1),l,1)
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
        end do
    end do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 2nd goup of direction cx<0 & cy>0 & cz>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    do l=1,Nc8
        do i=1,Nstencil2
            !k=di1(i)
            ii = di2(1,i)
            jj = di2(2,i)
            kk = di2(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch fom reflected in x-direction (default) to in y-direction
            ! only fo group of velocity overlapped by two reflected x&y-direction
            if (image(ii,jj-1,kk)==wallWN) then
                f(k-Nxtotal,l,2)=f(k-Nxtotal, l,4)
            elseif ((image(ii,jj-1,kk)==wallWNF) &
               .o. (image(ii,jj-1,kk)==wallWNB)) then
                f(k-Nxtotal,l,2)=fw(which_coner(ii,jj-1,kk),l,2)
            endif
            ! Switch fom reflected in x/y-direction to in z-direction
            ! only fo group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk-1)==wallNF) then
                f(k-Nxytotal,l,2)=f(k-Nxytotal,l,7)
            elseif (image(ii,jj,kk-1)==wallWF) then
                f(k-Nxytotal,l,2)=f(k-Nxytotal,l,5)
            elseif ((image(ii,jj,kk-1)==wallENF) &
               .o. (image(ii,jj,kk-1)==wallWSF)) then
                f(k-Nxytotal,l,2)=fw(which_coner(ii,jj,kk-1),l,2)
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
        end do
    end do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 3d group of direction cx<0 & cy<0 & cz>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    do l=1,Nc8
        do i=1,Nstencil3
            ii = di3(1,i)
            jj = di3(2,i)
            kk = di3(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch fom reflected in x-direction (default) to in y-direction
            ! only fo group of velocity overlapped by two reflected direction
            if (image(ii,jj+1,kk)==wallWS) then
                f(k+Nxtotal,l,3)=f(k+Nxtotal,l,1)
            elseif ((image(ii,jj+1,kk)==wallWSF) &
               .o. (image(ii,jj+1,kk)==wallWSB)) then
                f(k+Nxtotal,l,3)=fw(which_coner(ii,jj+1,kk),l,3)
            end if
            ! Switch fom reflected in x/y-direction to in z-direction
            ! only fo group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk-1)==wallSF) then
                f(k-Nxytotal,l,3)=f(k-Nxytotal,l,6) 
            elseif (image(ii,jj,kk-1)==wallWF) then
                f(k-Nxytotal,l,3)=f(k-Nxytotal,l,8) 
            elseif ((image(ii,jj,kk-1)==wallWNF) &
               .o. (image(ii,jj,kk-1)==wallESF)) then
                f(k-Nxytotal,l,3)=fw(which_coner(ii,jj,kk-1),l,3)
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
        end do
    end do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 4th goup of direction cx>0 & cy<0 & cz>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    do l=1,Nc8
        do i=1,Nstencil4
            ii = di4(1,i)
            jj = di4(2,i)
            kk = di4(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch fom reflected in x-direction (default) to in y-direction
            ! only fo group of velocity overlapped by two reflected direction
            if (image(ii,jj+1,kk)==wallES) then
                f(k+Nxtotal,l,4)=f(k+Nxtotal,l,2)
            elseif ((image(ii,jj+1,kk)==wallESF) &
               .o. (image(ii,jj+1,kk)==wallESB)) then
                f(k+Nxtotal,l,4)=fw(which_coner(ii,jj+1,kk),l,4)
            endif
            ! Switch fom reflected in x/y-direction to in z-direction
            ! only fo group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk-1)==wallSF) then
                f(k-Nxytotal,l,4)=f(k-Nxytotal,l,5)
            elseif (image(ii,jj,kk-1)==wallEF) then
                f(k-Nxytotal,l,4)=f(k-Nxytotal,l,7)
            elseif ((image(ii,jj,kk-1)==wallENF) &
                .o.(image(ii,jj,kk-1)==wallWSF)) then
                f(k-Nxytotal,l,4)=fw(which_coner(ii,jj,kk-1),l,4)
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
        end do
    end do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 5th goup of direction cx>0 & cy>0 & cz<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    do l=1,Nc8
        do i=1,Nstencil5
            ii = di5(1,i)
            jj = di5(2,i)
            kk = di5(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch fom reflected in x-direction (default) to in y-direction
            ! only fo group of velocity overlapped by two reflected direction x-y
            if (image(ii,jj-1,kk)==wallEN) then
                f(k-Nxtotal,l,5)=f(k-Nxtotal,l,7)
            elseif ((image(ii,jj-1,kk)==wallENF) &
               .o. (image(ii,jj-1,kk)==wallENB)) then
                f(k-Nxtotal,l,5)=fw(which_coner(ii,jj-1,kk),l,5)
            endif
            ! Switch fom reflected in x/y-direction to in z-direction
            ! only fo group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk+1)==wallNB) then
                f(k+Nxytotal,l,5)=f(k+Nxytotal,l,4)
            elseif (image(ii,jj,kk+1)==wallEB) then
                f(k+Nxytotal,l,5)=f(k+Nxytotal,l,2)
            elseif ((image(ii,jj,kk+1)==wallESB) &
               .o. (image(ii,jj,kk+1)==wallWNB)) then
                f(k+Nxytotal,l,5)=fw(which_coner(ii,jj,kk+1),l,5)
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
        end do
    end do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 6th goup of direction cx<0 & cy>0 & cz<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    do l=1,Nc8
        do i=1,Nstencil6
            ii = di6(1,i)
            jj = di6(2,i)
            kk = di6(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch fom reflected in x-direction (default) to in y-direction
            ! only fo group of velocity overlapped by two reflected x&y-direction
            if (image(ii,jj-1,kk)==wallWN) then
                f(k-Nxtotal,l,6)=f(k-Nxtotal,l,8)
            elseif ((image(ii,jj-1,kk)==wallWNF) &
                .o.(image(ii,jj-1,kk)==wallWNB)) then
                f(k-Nxtotal,l,6)=fw(which_coner(ii,jj-1,kk),l,6)
            endif
            ! Switch fom reflected in x/y-direction to in z-direction
            ! only fo group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk+1)==wallNB) then
                f(k+Nxytotal,l,6)=f(k+Nxytotal,l,3)
            elseif (image(ii,jj,kk+1)==wallWB) then
                f(k+Nxytotal,l,6)=f(k+Nxytotal,l,1)
            elseif ((image(ii,jj,kk+1)==wallENB) &
               .o. (image(ii,jj,kk+1)==wallWSB)) then
                f(k+Nxytotal,l,6)=fw(which_coner(ii,jj,kk+1),l,6)
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
        end do
    end do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 7th goup of direction cx<0 & cy<0 & cz<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    do l=1,Nc8
        do i=1,Nstencil7
            ii = di7(1,i)
            jj = di7(2,i)
            kk = di7(3,i) 
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch fom reflected in x-direction (default) to in y-direction
            ! only fo group of velocity overlapped by two reflected direction
            if (image(ii,jj+1,kk)==wallWS) then
                f(k+Nxtotal,l,7)=f(k+Nxtotal,l,8)
            elseif ((image(ii,jj+1,kk)==wallWSF) &
               .o. (image(ii,jj+1,kk)==wallWSB)) then
                f(k+Nxtotal,l,7)=fw(which_coner(ii,jj+1,kk),l,7)
            endif
            ! Switch fom reflected in x/y-direction to in z-direction
            ! only fo group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk+1)==wallSB) then
                f(k+Nxytotal,l,7)=f(k+Nxytotal,l,2)
            elseif (image(ii,jj,kk+1)==wallWB) then
                f(k+Nxytotal,l,7)=f(k+Nxytotal,l,4)
            elseif ((image(ii,jj,kk+1)==wallWNB) &
               .o. (image(ii,jj,kk+1)==wallESB)) then
                f(k+Nxytotal,l,7)=fw(which_coner(ii,jj,kk+1),l,7)
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
        end do
    end do
!$OMP END DO NOWAIT

!------------------------------------------------------------------------
!           In the 8th goup of direction cx>0 & cy<0 & cz<0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    do l=1,Nc8
        do i=1,Nstencil8
            ii = di8(1,i)
            jj = di8(2,i)
            kk = di8(3,i)
            k = (ii-xlg+1) + (jj-ylg)*Nxtotal + (kk-zlg)*Nxytotal
            ! Switch fom reflected in x-direction (default) to in y-direction
            ! only fo group of velocity overlapped by two reflected direction
            if (image(ii,jj+1,kk)==wallES) then
                f(k+Nxtotal,l,8)=f(k+Nxtotal,l,6)
            elseif ((image(ii,jj+1,kk)==wallESF) &
               .o. (image(ii,jj+1,kk)==wallESB)) then
                f(k+Nxtotal,l,8)=fw(which_coner(ii,jj+1,kk),l,8)
            endif
            ! Switch fom reflected in x/y-direction to in z-direction
            ! only fo group of velocity overlapped by two/three reflected x/y&z-direction
            if (image(ii,jj,kk+1)==wallSB) then
                f(k+Nxytotal,l,8)=f(k+Nxytotal,l,1)
            elseif (image(ii,jj,kk+1)==wallEB) then
                f(k+Nxytotal,l,8)=f(k+Nxytotal,l,3)
            elseif ((image(ii,jj,kk+1)==wallENB) &
               .o. (image(ii,jj,kk+1)==wallWSB)) then
                f(k+Nxytotal,l,8)=fw(which_coner(ii,jj,kk+1),l,8)
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
        end do
    end do
!$OMP END DO 

!$OMP SINGLE
        ! Wait until send and ecv done
        call MPI_WAITALL(4, MPI_REQ_X, MPI_STAT, MPI_ERR)
        call MPI_WAITALL(4, MPI_REQ_Y, MPI_STAT, MPI_ERR)
        call MPI_WAITALL(4, MPI_REQ_Z, MPI_STAT, MPI_ERR)
!$OMP END SINGLE 

!$OMP DO COLLAPSE(3)
    ! pack/unpack X di buffer
    do k = 1,Nztotal
        do j = 1, Nytotal
            do i = 1, ghostLayes
                ! pack oder, 
                !fo dir: (2,3,6,7)
                localidMsnd = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i+ghostLayes
                localidMcv = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i
                !fo dir: (1,4,5,8)
                localidPsnd = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i+Nxsub
                localidPcv = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i+ghostLayers+Nxsub
        
                packid = (k-1)*Nytotal*ghostLayes*Nc &
                       + (j-1)*ghostLayes*Nc &
                       + (i-1)*Nc
                do l = 1,8
                    f_west_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidMsnd,:,l)
                    f(localidMcv,:,l) = f_west_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                    f_east_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidPsnd,:,l)
                    f(localidPcv,:,l) = f_east_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                enddo !l
            enddo !i
        enddo !j
    enddo !k
!$OMP END DO NOWAIT

!$OMP DO COLLAPSE(3)
    ! pack/unpack Y di buffer
    do k = 1, Nztotal
        do j = 1, ghostLayes
            do i = 1, Nxtotal
                !fo dir: (3,4,7,8)
                localidMsnd = (k-1)*Nxtotal*Nytotal + (j+ghostLayes-1)*Nxtotal + i
                localidMcv = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i
                !fo dir: (2,1,6,5)
                localidPsnd = (k-1)*Nxtotal*Nytotal + (j+Nysub-1)*Nxtotal + i
                localidPcv = (k-1)*Nxtotal*Nytotal + (j+Nysub+ghostLayers-1)*Nxtotal + i
                packid = (k-1)*ghostLayes*Nxtotal*Nc &
                       + (j-1)*Nxtotal*Nc &
                       + (i-1)*Nc
                do l = 1,8
                    f_suth_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidMsnd,:,l)
                    f(localidMcv,:,l) = f_suth_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                    f_noth_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidPsnd,:,l)
                    f(localidPcv,:,l) = f_noth_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                enddo !l
            enddo !i
        enddo !j
    enddo !k
!$OMP END DO NOWAIT


!$OMP DO COLLAPSE(3)
    ! pack/unpack Z di buffer
    do k = 1, ghostLayes
        do j = 1, Nytotal
            do i = 1, Nxtotal
                ! pack oder, 
                !fo dir: (5,6,7,8)
                localidMsnd = (k+ghostLayes-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                localidMcv = (k-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                !fo dir: (1,2,3,4)
                localidPsnd = (k+Nzsub-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                localidPcv = (k+Nzsub+ghostLayers-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                packid = (k-1)*Nytotal*Nxtotal*Nc &
                       + (j-1)*Nxtotal*Nc &
                       + (i-1)*Nc
                do l = 1,8
                    f_back_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidMsnd,:,l)
                    f(localidMcv,:,l) = f_back_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                    f_fnt_snd(packid+1+(l-1)*Nc8:packid+l*Nc8) = f(localidPsnd,:,l)
                    f(localidPcv,:,l) = f_frnt_rcv(packid+1+(l-1)*Nc8:packid+l*Nc8)
                enddo !l                
            enddo !i
        enddo !j
    enddo !k
!$OMP END DO NOWAIT


!------------------------------------------------------------------------
! pack/unpack 3-fold coners at west/east
!------------------------------------------------------------------------
!$OMP SECTIONS
!$OMP SECTION
!!$OMP DO
    ! pack 3-fold coners send buffer at west
    do m = 1, westN3coner_snd
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            f_west_snd(locB+1:locB+Nc8) = fw(map3CoWsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! pack 3-fold coners send buffer at east
    do m = 1, eastN3coner_snd
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            f_east_snd(locB+1:locB+Nc8) = fw(map3CoEsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold coners rcv buffer at west
    do m = 1, westN3coner_rcv
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CoWrcv(m),:,l) = f_west_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold coners rcv buffer at east
    do m = 1, eastN3coner_rcv
        do l = 1,8
            locB = xfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CoErcv(m),:,l) = f_east_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO NOWAIT

!------------------------------------------------------------------------
! pack/unpack 3-fold coners at noth/suth
!------------------------------------------------------------------------
!$OMP SECTION
!!$OMP DO
    ! pack 3-fold coners send buffer at suth
    do m = 1, suthN3coner_snd
        do l = 1,8
            locB = yfsize+(m-1)*Nc+(l-1)*Nc8
            f_suth_snd(locB+1:locB+Nc8) = fw(map3CoSsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO

!$OMP SECTION
!!$OMP DO
    ! pack 3-fold coners send buffer at noth
    do m = 1, nothN3coner_snd
        do l = 1,8
            locB = yfsize+(m-1)*Nc+(l-1)*Nc8
            f_noth_snd(locB+1:locB+Nc8) = fw(map3CoNsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold coners rcv buffer at suth
    do m = 1, suthN3coner_rcv
        do l = 1,8
            locB = yfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CoSrcv(m),:,l) = f_suth_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold coners rcv buffer at noth
    do m = 1, nothN3coner_rcv
        do l = 1,8
            locB = yfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CoNrcv(m),:,l) = f_noth_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO

!------------------------------------------------------------------------
! pack/unpack 3-fold coners at back/frnt
!------------------------------------------------------------------------
!$OMP SECTION
!!$OMP DO
    ! pack 3-fold coners send buffer at back
    do m = 1, backN3coner_snd
        do l = 1,8
            locB = zfsize+(m-1)*Nc+(l-1)*Nc8
            f_back_snd(locB+1:locB+Nc8) = fw(map3CoBsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! pack 3-fold coners send buffer at frnt
    do m = 1, fntN3corner_snd
        do l = 1,8
            locB = zfsize+(m-1)*Nc+(l-1)*Nc8
            f_fnt_snd(locB+1:locB+Nc8) = fw(map3CorFsnd(m),:,l)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold coners rcv buffer at back
    do m = 1, backN3coner_rcv
        do l = 1,8
            locB = zfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CoBrcv(m),:,l) = f_back_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO NOWAIT

!$OMP SECTION
!!$OMP DO
    ! unpack 3-fold coners rcv buffer at frnt
    do m = 1, fntN3corner_rcv
        do l = 1,8
            locB = zfsize+(m-1)*Nc+(l-1)*Nc8
            fw(map3CoFrcv(m),:,l) = f_frnt_rcv(locB+1:locB+Nc8)
        end do
    end do
!!$OMP ENDDO NOWAIT
!$OMP END SECTIONS

!=======================================================================
!     Bounday condition on the flat wall
!=======================================================================
!$OMP DO SCHEDULE(STATIC)
    do i=1,nWall
        k=vecWall(i)
        kk = k/Nxytotal + zlg ! to be checked
        jj = (k-(kk-zlg)*Nxytotal)/Nxtotal + ylg
        ii = k-(kk-zlg)*Nxytotal-(jj-ylg)*Nxtotal + xlg -1      
        j=which_coner(ii,jj,kk)
        RhoWall=0.d0
        RhoWall2=0.d0
        RhoWall3=0.d0
        select case (image(ii,jj,kk))
            case (wallE)
                do l=1,Nc8
                    f(k,l,2)=extCoef(i,1)*f(k+1,l,2)+extCoef(i,2)*f(k+2,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k+1,l,3)+extCoef(i,2)*f(k+2,l,3)
                    f(k,l,6)=extCoef(i,1)*f(k+1,l,6)+extCoef(i,2)*f(k+2,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+1,l,7)+extCoef(i,2)*f(k+2,l,7)
                    RhoWall=RhoWall+cx(l)*(f(k,l,2)+f(k,l,3)+f(k,l,6)+f(k,l,7))
                enddo
                RhoWall=RhoWall/DiffFlux
                do l=1,Nc8
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                enddo
            case (wallW)
                do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-1,l,1)+extCoef(i,2)*f(k-2,l,1)
                    f(k,l,4)=extCoef(i,1)*f(k-1,l,4)+extCoef(i,2)*f(k-2,l,4)
                    f(k,l,5)=extCoef(i,1)*f(k-1,l,5)+extCoef(i,2)*f(k-2,l,5)
                    f(k,l,8)=extCoef(i,1)*f(k-1,l,8)+extCoef(i,2)*f(k-2,l,8)
                    RhoWall=RhoWall+cx(l)*(f(k,l,1)+f(k,l,4)+f(k,l,5)+f(k,l,8))
                enddo
                RhoWall=RhoWall/DiffFlux
                do l=1,Nc8
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                enddo
            case (wallN)
                do l=1,Nc8
                    f(k,l,3)=extCoef(i,1)*f(k+Nxtotal,l,3)+extCoef(i,2)*f(k+2*Nxtotal,l,3)
                    f(k,l,4)=extCoef(i,1)*f(k+Nxtotal,l,4)+extCoef(i,2)*f(k+2*Nxtotal,l,4)
                    f(k,l,7)=extCoef(i,1)*f(k+Nxtotal,l,7)+extCoef(i,2)*f(k+2*Nxtotal,l,7)
                    f(k,l,8)=extCoef(i,1)*f(k+Nxtotal,l,8)+extCoef(i,2)*f(k+2*Nxtotal,l,8)
                    RhoWall=RhoWall+cy(l)*(f(k,l,3)+f(k,l,4)+f(k,l,7)+f(k,l,8))
                enddo
                RhoWall=RhoWall/DiffFlux
                do l=1,Nc8
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                enddo
            case (wallS)
                do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-Nxtotal,l,1)+extCoef(i,2)*f(k-2*Nxtotal,l,1)
                    f(k,l,2)=extCoef(i,1)*f(k-Nxtotal,l,2)+extCoef(i,2)*f(k-2*Nxtotal,l,2)
                    f(k,l,5)=extCoef(i,1)*f(k-Nxtotal,l,5)+extCoef(i,2)*f(k-2*Nxtotal,l,5)
                    f(k,l,6)=extCoef(i,1)*f(k-Nxtotal,l,6)+extCoef(i,2)*f(k-2*Nxtotal,l,6)
                    RhoWall=RhoWall+cy(l)*(f(k,l,1)+f(k,l,2)+f(k,l,5)+f(k,l,6))
                enddo
                RhoWall=RhoWall/DiffFlux
                do l=1,Nc8
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                enddo
            case (wallF)
                do l=1,Nc8
                    f(k,l,5)=extCoef(i,1)*f(k+Nxytotal,l,5)+extCoef(i,2)*f(k+2*Nxytotal,l,5)
                    f(k,l,6)=extCoef(i,1)*f(k+Nxytotal,l,6)+extCoef(i,2)*f(k+2*Nxytotal,l,6)
                    f(k,l,7)=extCoef(i,1)*f(k+Nxytotal,l,7)+extCoef(i,2)*f(k+2*Nxytotal,l,7)
                    f(k,l,8)=extCoef(i,1)*f(k+Nxytotal,l,8)+extCoef(i,2)*f(k+2*Nxytotal,l,8)
                    RhoWall=RhoWall+cz(l)*(f(k,l,5)+f(k,l,6)+f(k,l,7)+f(k,l,8))
                enddo
                RhoWall=RhoWall/DiffFlux
                do l=1,Nc8
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                enddo
            case (wallB)
                do l=1,Nc8
                    f(k,l,1)=extCoef(i,1)*f(k-Nxytotal,l,1)+extCoef(i,2)*f(k-2*Nxytotal,l,1)
                    f(k,l,2)=extCoef(i,1)*f(k-Nxytotal,l,2)+extCoef(i,2)*f(k-2*Nxytotal,l,2)
                    f(k,l,3)=extCoef(i,1)*f(k-Nxytotal,l,3)+extCoef(i,2)*f(k-2*Nxytotal,l,3)
                    f(k,l,4)=extCoef(i,1)*f(k-Nxytotal,l,4)+extCoef(i,2)*f(k-2*Nxytotal,l,4)
                    RhoWall=RhoWall+cz(l)*(f(k,l,1)+f(k,l,2)+f(k,l,3)+f(k,l,4))
                enddo
                RhoWall=RhoWall/DiffFlux
                do l=1,Nc8
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                enddo
!=======================================================================
!     Bounday condition on the 2-direction corner wall
!=======================================================================
!------------------------------------------------------------------------
!           No z-diection type wall
!------------------------------------------------------------------------
            case (wallEN)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite y
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    ! Stoe y to unused location
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                enddo

            case (wallWN)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8

                    ! Wite x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Wite y
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Stoe y to unused location
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                enddo

            case (wallES)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite y
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    ! Stoe y to unused location
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                enddo

            case (wallWS)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Wite y
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    ! Stoe y to unused location
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                enddo

!------------------------------------------------------------------------
!           No y-diection type wall
!------------------------------------------------------------------------
            case (wallEF)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite z
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    ! Stoe z to unsed location
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                enddo

            case (wallWF)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Wite z
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Stoe z to unsed location
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                enddo

            case (wallEB)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite z
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    ! Stoe z to unsed location
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                enddo

            case (wallWB)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Wite z
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    ! Stoe z
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                enddo
!------------------------------------------------------------------------
!           No x-diection type wall
!------------------------------------------------------------------------
            case (wallNF)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite y
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite z
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Stoe z to unsed location
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                enddo

            case (wallSF)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite y
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    ! Wite z
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    ! Stoe z to unsed location
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,8)
                enddo

            case (wallNB)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite y
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite z
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    ! Stoe z to unsed location
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                enddo

            case (wallSB)
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                do l=1,Nc8
                    ! Wite y
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    ! Wite z
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    ! Stoe z to unsed location
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fwZ(l,4)
                enddo
!=======================================================================
!     Bounday condition on the 3-direction corner wall
!=======================================================================
            case (wallENF) !diection1
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite y
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,3)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,7)
                    ! Wite z
                    f(k,l,3)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,7)

                    ! Stoe y
                    fw(j,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,4)
                    fw(j,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,8)
                    ! Stoe z
                    fw(j,l,2)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    fw(j,l,4)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Borow the un-used variable to store z for the opposite velocity
                    f(k,l,7)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,5)
                enddo
            case (wallWNF) !diection2
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Wite y
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,4)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,8)
                    ! Wite z
                    f(k,l,4)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,8)

                    ! Stoe y
                    fw(j,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,3)
                    fw(j,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,7)
                    ! Stoe z
                    fw(j,l,1)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    fw(j,l,3)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    ! Borow the un-used variable to store z for the opposite velocity
                    f(k,l,8)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,6)
                enddo
            case (wallWSF) !diection3
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Wite y
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,1)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,5)
                    ! Wite z
                    f(k,l,1)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,5)

                    ! Stoe y
                    fw(j,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,2)
                    fw(j,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,6)
                    ! Stoe z
                    fw(j,l,2)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,6)
                    fw(j,l,4)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,8)
                    ! Borow the un-used variable to store z for the opposite velocity
                    f(k,l,5)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,7)
                enddo
            case (wallESF) !diection4
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite y
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,2)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,6)
                    ! Wite z
                    f(k,l,2)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,6)

                    ! Stoe y
                    fw(j,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,1)
                    fw(j,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,5)
                    ! Stoe z
                    fw(j,l,1)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,5)
                    fw(j,l,3)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,7)
                    ! Borow the un-used variable to store z for the opposite velocity
                    f(k,l,6)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,8)
                enddo
            case (wallENB) !diection5
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite y
                    f(k,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,3)
                    f(k,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,7)
                    ! Wite z
                    f(k,l,7)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,3)

                    ! Stoe y
                    fw(j,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,4)
                    fw(j,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,8)
                    ! Stoe z
                    fw(j,l,6)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    fw(j,l,8)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    ! Borow the un-used variable to store z for the opposite velocity
                    f(k,l,3)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,1)
                enddo
            case (wallWNB) !diection6
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Wite y
                    f(k,l,1)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,4)
                    f(k,l,5)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,8)
                    ! Wite z
                    f(k,l,8)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,4)

                    ! Stoe y
                    fw(j,l,2)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,3)
                    fw(j,l,6)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,7)
                    ! Stoe z
                    fw(j,l,5)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    fw(j,l,7)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    ! Borow the un-used variable to store z for the opposite velocity
                    f(k,l,4)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,2)
                enddo
            case (wallWSB) !diection7
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,2)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,1)
                    f(k,l,3)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,4)
                    f(k,l,6)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,5)
                    f(k,l,7)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,8)
                    ! Wite y
                    f(k,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,1)
                    f(k,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,5)
                    ! Wite z
                    f(k,l,5)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,1)

                    ! Stoe y
                    fw(j,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,2)
                    fw(j,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,6)
                    ! Stoe z
                    fw(j,l,6)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,2)
                    fw(j,l,8)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,4)
                    ! Borow the un-used variable to store z for the opposite velocity
                    f(k,l,1)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,3)
                enddo
            case (wallESB) !diection8
                do l=1,Nc8
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
                enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                do l=1,Nc8
                    ! Wite x
                    f(k,l,1)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,2)
                    f(k,l,4)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,3)
                    f(k,l,5)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,6)
                    f(k,l,8)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f(k,l,7)
                    ! Wite y
                    f(k,l,3)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,2)
                    f(k,l,7)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,6)
                    ! Wite z
                    f(k,l,6)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,2)

                    ! Stoe y
                    fw(j,l,4)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,1)
                    fw(j,l,8)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*fw(j,l,5)
                    ! Stoe z
                    fw(j,l,5)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,1)
                    fw(j,l,7)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,3)
                    ! Borow the un-used variable to store z for the opposite velocity
                    f(k,l,2)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*fwZ(l,4)
                enddo
        end select
    end do
!$OMP END DO   


!--------------------------------------------------------
!> inlet/outlet
!--------------------------------------------------------
    if(xl==xmin) then ! inlet block (west most pocessor)
!$OMP DO
        do k=zlg,zug
            do j=ylg,yug
                i = xl
                l = (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                !inlet
                f(l,:,1)=f(l-1,:,1) + w(:)*PessDrop
                f(l,:,4)=f(l-1,:,4) + w(:)*PessDrop
                f(l,:,5)=f(l-1,:,5) + w(:)*PessDrop
                f(l,:,8)=f(l-1,:,8) + w(:)*PessDrop
            end do
        end do
!$OMP END DO
    endif

    if(xu==xmax) then ! outlet block (east most pocessor)
!$OMP DO
        do k=zlg,zug
            do j=ylg,yug
                i = xu
                l = (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                !outlet
                f(l,:,2)=f(l+1,:,2) - w(:)*PessDrop
                f(l,:,3)=f(l+1,:,3) - w(:)*PessDrop
                f(l,:,6)=f(l+1,:,6) - w(:)*PessDrop
                f(l,:,7)=f(l+1,:,7) - w(:)*PessDrop
            enddo
        end do
!$OMP END DO            
    endif

!----------------------------------------------------
!> Symmetic BC
!----------------------------------------------------
    if(yl==ymin) then ! south sym BC
!$OMP DO
        do k=zl,zu
            do i=xl,xu
                l = (k-zlg)*Nxytotal + ghostLayes*Nxtotal + i-xlg+1
                f(l,:,2)=f(l,:,3)
                f(l,:,1)=f(l,:,4)
                f(l,:,6)=f(l,:,7)
                f(l,:,5)=f(l,:,8)
            enddo
        end do
!$OMP END DO 
    endif
    if(yu==ymax) then ! noth sym BC
!$OMP DO
        do k=zl,zu
            do i=xl,xu
                l = (k-zlg)*Nxytotal + (ghostLayes+Nysub-1)*Nxtotal + i-xlg+1
                f(l,:,3)=f(l,:,2)
                f(l,:,4)=f(l,:,1)
                f(l,:,7)=f(l,:,6)
                f(l,:,8)=f(l,:,5)
            end do
        end do
!$OMP END DO
    endif
    if(zl==zmin) then ! back sym BC
!$OMP DO
        do j=yl,yu
            do i=xl,xu
                l = ghostLayes*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                f(l,:,1)=f(l,:,5)
                f(l,:,2)=f(l,:,6)
                f(l,:,3)=f(l,:,7)
                f(l,:,4)=f(l,:,8)
            enddo
        end do
!$OMP END DO 
    endif
    if(zu==zmax) then ! font sym BC
!$OMP DO
        do j=yl,yu
            do i=xl,xu
                l = (ghostLayes+Nzsub-1)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                f(l,:,5)=f(l,:,1)
                f(l,:,6)=f(l,:,2)
                f(l,:,7)=f(l,:,3)
                f(l,:,8)=f(l,:,4)
            end do
        end do
!$OMP END DO
    endif

!----------------------------------------------------
!> Update Maco
!----------------------------------------------------
!$OMP DO SCHEDULE(STATIC) 
    !do k=1,Ntotal
    do i=1,Nfluid
        k=mapF(i)
        Rho(k)=0.d0
        Ux(k)=0.d0
        Uy(k)=0.d0
        Uz(k)=0.d0
        do l=1,Nc8
            Rho(k)=Rho(k)+f(k,l,1)+f(k,l,2)+f(k,l,3)+f(k,l,4)+f(k,l,5)+f(k,l,6)+f(k,l,7)+f(k,l,8)
            Ux(k)=Ux(k)+cx(l)*(f(k,l,1)-f(k,l,2)-f(k,l,3)+f(k,l,4)+f(k,l,5)-f(k,l,6)-f(k,l,7)+f(k,l,8))
            Uy(k)=Uy(k)+cy(l)*(f(k,l,1)+f(k,l,2)-f(k,l,3)-f(k,l,4)+f(k,l,5)+f(k,l,6)-f(k,l,7)-f(k,l,8))
            Uz(k)=Uz(k)+cz(l)*(f(k,l,1)+f(k,l,2)+f(k,l,3)+f(k,l,4)-f(k,l,5)-f(k,l,6)-f(k,l,7)-f(k,l,8))
        end do
    end do
!$OMP END DO 
!$OMP END PARALLEL  
    end suboutine iterate
end module solve
