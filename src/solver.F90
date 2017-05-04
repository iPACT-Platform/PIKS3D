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

        integer :: i, j, k, l, ii, jj, kk
        integer :: localidMsnd, localidPsnd, localidMrcv, localidPrcv, packid
        INTEGER :: MPI_ERR
        INTEGER :: MPI_REQ_X(4), MPI_REQ_Y(4), MPI_REQ_Z(4)
        INTEGER :: MPI_STAT(MPI_STATUS_SIZE,6)
        integer :: xsize, ysize, zsize
        double precision :: fEq, RhoWall, RhoWall2, RhoWall3
        double precision, dimension(Nc8) :: f1wZ,f2wZ,f3wZ,f4wZ,f5wZ,f6wZ,f7wZ,f8wZ

        ! buffer size
        xsize = Nytotal*Nztotal*Nc/2*ghostLayers
        ysize = Nxtotal*Nztotal*Nc/2*ghostLayers
        zsize = Nxtotal*Nytotal*Nc/2*ghostLayers

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(l, i, j, k, ii, jj, kk, fEq, RhoWall, RhoWall2, RhoWall3)&
!$OMP PRIVATE(f1wZ,f2wZ,f3wZ,f4wZ,f5wZ,f6wZ,f7wZ,f8wZ)&
!$OMP PRIVATE(localidPsnd, localidPrcv, localidMsnd, localidMrcv, packid)


!$OMP SINGLE
        ! Start Recieving
        CALL MPI_IRECV( f_east_rcv, xsize, MPI_DOUBLE_PRECISION, east, TAG1, &
                        MPI_COMM_VGRID, MPI_REQ_X(1), MPI_ERR )
        CALL MPI_IRECV( f_west_rcv, xsize, MPI_DOUBLE_PRECISION, west, TAG2, &
                        MPI_COMM_VGRID, MPI_REQ_X(2), MPI_ERR )
        CALL MPI_IRECV( f_noth_rcv, ysize, MPI_DOUBLE_PRECISION, noth, TAG3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(1), MPI_ERR )
        CALL MPI_IRECV( f_suth_rcv, ysize, MPI_DOUBLE_PRECISION, suth, TAG4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(2), MPI_ERR )     
        CALL MPI_IRECV( f_frnt_rcv, zsize, MPI_DOUBLE_PRECISION, frnt, TAG5, &
                        MPI_COMM_VGRID, MPI_REQ_Z(1), MPI_ERR )
        CALL MPI_IRECV( f_back_rcv, zsize, MPI_DOUBLE_PRECISION, back, TAG6, &
                        MPI_COMM_VGRID, MPI_REQ_Z(2), MPI_ERR )   
        ! Start Sending
        CALL MPI_ISEND( f_west_snd, xsize, MPI_DOUBLE_PRECISION, west, TAG1, &
                        MPI_COMM_VGRID, MPI_REQ_X(3), MPI_ERR )
        CALL MPI_ISEND( f_east_snd, xsize, MPI_DOUBLE_PRECISION, east, TAG2, &
                        MPI_COMM_VGRID, MPI_REQ_X(4), MPI_ERR )
        CALL MPI_ISEND( f_suth_snd, ysize, MPI_DOUBLE_PRECISION, suth, TAG3, &
                        MPI_COMM_VGRID, MPI_REQ_Y(3), MPI_ERR )
        CALL MPI_ISEND( f_noth_snd, ysize, MPI_DOUBLE_PRECISION, noth, TAG4, &
                        MPI_COMM_VGRID, MPI_REQ_Y(4), MPI_ERR ) 
        CALL MPI_ISEND( f_back_snd, zsize, MPI_DOUBLE_PRECISION, back, TAG5, &
                        MPI_COMM_VGRID, MPI_REQ_Z(3), MPI_ERR )
        CALL MPI_ISEND( f_frnt_snd, zsize, MPI_DOUBLE_PRECISION, frnt, TAG6, &
                        MPI_COMM_VGRID, MPI_REQ_Z(4), MPI_ERR ) 

!$OMP END SINGLE NOWAIT

!------------------------------------------------------------------------
!           In the 1st group of direction cx<0 & cy>0 & cz>0
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
            If ((image(ii,jj-1,kk)==wallEN).OR.(image(ii,jj-1,kk)==wallENF) &
            .OR.(image(ii,jj-1,kk)==wallENB)) then
                f1(k-Nxtotal,l)=f1w(which_corner(ii,jj-1,kk),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(ii,jj,kk-1)==wallNF) .OR.(image(ii,jj,kk-1)==wallEF) &
            .OR.(image(ii,jj,kk-1)==wallWNF).OR.(image(ii,jj,kk-1)==wallESF)) then
                f1(k-Nxytotal,l)=f1w(which_corner(ii,jj,kk-1),l)
            else if (image(ii,jj,kk-1)==wallENF) then
                f1(k-Nxytotal,l)=f7(k-Nxytotal,l)
            End if

            fEq=w(l)*(Rho(k)+1.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)+cz(l)*Uz(k)))
            f1(k,l)=(mu*(fEq-0.5d0*f1(k,l)) &
            &        + cx(l)*coef1(i,2)*f1(k-1,l) &
            &        + cx(l)*coef1(i,3)*f1(k-2,l) &
            &        + cy(l)*coef1(i,5)*f1(k-Nxtotal,l) &
            &        + cy(l)*coef1(i,6)*f1(k-2*Nxtotal,l) &
            &        + cz(l)*coef1(i,8)*f1(k-Nxytotal,l) &
            &        + cz(l)*coef1(i,9)*f1(k-2*Nxytotal,l) &
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
            If ((image(ii,jj-1,kk)==wallWN).OR.(image(ii,jj-1,kk)==wallWNF) &
            .OR.(image(ii,jj-1,kk)==wallWNB)) then
                f2(k-Nxtotal,l)=f2w(which_corner(ii,jj-1,kk),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(ii,jj,kk-1)==wallNF) .OR.(image(ii,jj,kk-1)==wallWF) &
            .OR.(image(ii,jj,kk-1)==wallENF).OR.(image(ii,jj,kk-1)==wallWSF)) then
                f2(k-Nxytotal,l)=f2w(which_corner(ii,jj,kk-1),l)
            else if (image(ii,jj,kk-1)==wallWNF) then
                f2(k-Nxytotal,l)=f8(k-Nxytotal,l)
            End if

            fEq=w(l)*(Rho(k)+1.d0*(-cx(l)*Ux(k)+cy(l)*Uy(k)+cz(l)*Uz(k)))
            f2(k,l)=(mu*(fEq-0.5d0*f2(k,l)) &
            &        - cx(l)*coef2(i,2)*f2(k+1,l) &
            &        - cx(l)*coef2(i,3)*f2(k+2,l) &
            &        + cy(l)*coef2(i,5)*f2(k-Nxtotal,l) &
            &        + cy(l)*coef2(i,6)*f2(k-2*Nxtotal,l) &
            &        + cz(l)*coef2(i,8)*f2(k-Nxytotal,l) &
            &        + cz(l)*coef2(i,9)*f2(k-2*Nxytotal,l) &
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
            If ((image(ii,jj+1,kk)==wallWS).OR.(image(ii,jj+1,kk)==wallWSF) &
            .OR.(image(ii,jj+1,kk)==wallWSB)) then
                f3(k+Nxtotal,l)=f3w(which_corner(ii,jj+1,kk),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(ii,jj,kk-1)==wallSF) .OR.(image(ii,jj,kk-1)==wallWF) &
            .OR.(image(ii,jj,kk-1)==wallWNF).OR.(image(ii,jj,kk-1)==wallESF)) then
                f3(k-Nxytotal,l)=f3w(which_corner(ii,jj,kk-1),l)
            else if (image(ii,jj,kk-1)==wallWSF) then
                f3(k-Nxytotal,l)=f5(k-Nxytotal,l)
            End if

            fEq=w(l)*(Rho(k)+1.d0*(-cx(l)*Ux(k)-cy(l)*Uy(k)+cz(l)*Uz(k)))
            f3(k,l)=(mu*(fEq-0.5d0*f3(k,l)) &
            &        - cx(l)*coef3(i,2)*f3(k+1,l) &
            &        - cx(l)*coef3(i,3)*f3(k+2,l) &
            &        - cy(l)*coef3(i,5)*f3(k+Nxtotal,l) &
            &        - cy(l)*coef3(i,6)*f3(k+2*Nxtotal,l) &
            &        + cz(l)*coef3(i,8)*f3(k-Nxytotal,l) &
            &        + cz(l)*coef3(i,9)*f3(k-2*Nxytotal,l) &
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
            If ((image(ii,jj+1,kk)==wallES).OR.(image(ii,jj+1,kk)==wallESF) &
            .OR.(image(ii,jj+1,kk)==wallESB)) then
                f4(k+Nxtotal,l)=f4w(which_corner(ii,jj+1,kk),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(ii,jj,kk-1)==wallSF) .OR.(image(ii,jj,kk-1)==wallEF) &
            .OR.(image(ii,jj,kk-1)==wallENF).OR.(image(ii,jj,kk-1)==wallWSF)) then
                f4(k-Nxytotal,l)=f4w(which_corner(ii,jj,kk-1),l)
            else if (image(ii,jj,kk-1)==wallESF) then
                f4(k-Nxytotal,l)=f6(k-Nxytotal,l)
            End if

            fEq=w(l)*(Rho(k)+1.d0*(cx(l)*Ux(k)-cy(l)*Uy(k)+cz(l)*Uz(k)))
            f4(k,l)=(mu*(fEq-0.5d0*f4(k,l)) &
            &        + cx(l)*coef4(i,2)*f4(k-1,l) &
            &        + cx(l)*coef4(i,3)*f4(k-2,l) &
            &        - cy(l)*coef4(i,5)*f4(k+Nxtotal,l) &
            &        - cy(l)*coef4(i,6)*f4(k+2*Nxtotal,l) &
            &        + cz(l)*coef4(i,8)*f4(k-Nxytotal,l) &
            &        + cz(l)*coef4(i,9)*f4(k-2*Nxytotal,l) &
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
            If ((image(ii,jj-1,kk)==wallEN).OR.(image(ii,jj-1,kk)==wallENF) &
            .OR.(image(ii,jj-1,kk)==wallENB)) then
                f5(k-Nxtotal,l)=f5w(which_corner(ii,jj-1,kk),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(ii,jj,kk+1)==wallNB) .OR.(image(ii,jj,kk+1)==wallEB) &
            .OR.(image(ii,jj,kk+1)==wallESB).OR.(image(ii,jj,kk+1)==wallWNB)) then
                f5(k+Nxytotal,l)=f5w(which_corner(ii,jj,kk+1),l)
            else if (image(ii,jj,kk+1)==wallENB) then
                f5(k+Nxytotal,l)=f3(k+Nxytotal,l)
            End if

            fEq=w(l)*(Rho(k)+1.d0*(cx(l)*Ux(k)+cy(l)*Uy(k)-cz(l)*Uz(k)))
            f5(k,l)=(mu*(fEq-0.5d0*f5(k,l)) &
            &        + cx(l)*coef5(i,2)*f5(k-1,l) &
            &        + cx(l)*coef5(i,3)*f5(k-2,l) &
            &        + cy(l)*coef5(i,5)*f5(k-Nxtotal,l) &
            &        + cy(l)*coef5(i,6)*f5(k-2*Nxtotal,l) &
            &        - cz(l)*coef5(i,8)*f5(k+Nxytotal,l) &
            &        - cz(l)*coef5(i,9)*f5(k+2*Nxytotal,l) &
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
            If ((image(ii,jj-1,kk)==wallWN).OR.(image(ii,jj-1,kk)==wallWNF) &
            .OR.(image(ii,jj-1,kk)==wallWNB)) then
                f6(k-Nxtotal,l)=f6w(which_corner(ii,jj-1,kk),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(ii,jj,kk+1)==wallNB) .OR.(image(ii,jj,kk+1)==wallWB) &
            .OR.(image(ii,jj,kk+1)==wallENB).OR.(image(ii,jj,kk+1)==wallWSB)) then
                f6(k+Nxytotal,l)=f6w(which_corner(ii,jj,kk+1),l)
            else if (image(ii,jj,kk+1)==wallWNB) then
                f6(k+Nxytotal,l)=f4(k+Nxytotal,l)
            End if

            fEq=w(l)*(Rho(k)+1.d0*(-cx(l)*Ux(k)+cy(l)*Uy(k)-cz(l)*Uz(k)))
            f6(k,l)=(mu*(fEq-0.5d0*f6(k,l)) &
            &        - cx(l)*coef6(i,2)*f6(k+1,l) &
            &        - cx(l)*coef6(i,3)*f6(k+2,l) &
            &        + cy(l)*coef6(i,5)*f6(k-Nxtotal,l) &
            &        + cy(l)*coef6(i,6)*f6(k-2*Nxtotal,l) &
            &        - cz(l)*coef6(i,8)*f6(k+Nxytotal,l) &
            &        - cz(l)*coef6(i,9)*f6(k+2*Nxytotal,l) &
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
            If ((image(ii,jj+1,kk)==wallWS).OR.(image(ii,jj+1,kk)==wallWSF) &
            .OR.(image(ii,jj+1,kk)==wallWSB)) then
                f7(k+Nxtotal,l)=f7w(which_corner(ii,jj+1,kk),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(ii,jj,kk+1)==wallSB) .OR.(image(ii,jj,kk+1)==wallWB) &
            .OR.(image(ii,jj,kk+1)==wallWNB).OR.(image(ii,jj,kk+1)==wallESB)) then
                f7(k+Nxytotal,l)=f7w(which_corner(ii,jj,kk+1),l)
            else if (image(ii,jj,kk+1)==wallWSB) then
                f7(k+Nxytotal,l)=f1(k+Nxytotal,l)
            End if

            fEq=w(l)*(Rho(k)+1.d0*(-cx(l)*Ux(k)-cy(l)*Uy(k)-cz(l)*Uz(k)))
            f7(k,l)=(mu*(fEq-0.5d0*f7(k,l)) &
            &        - cx(l)*coef7(i,2)*f7(k+1,l) &
            &        - cx(l)*coef7(i,3)*f7(k+2,l) &
            &        - cy(l)*coef7(i,5)*f7(k+Nxtotal,l) &
            &        - cy(l)*coef7(i,6)*f7(k+2*Nxtotal,l) &
            &        - cz(l)*coef7(i,8)*f7(k+Nxytotal,l) &
            &        - cz(l)*coef7(i,9)*f7(k+2*Nxytotal,l) &
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
            If ((image(ii,jj+1,kk)==wallES).OR.(image(ii,jj+1,kk)==wallESF) &
            .OR.(image(ii,jj+1,kk)==wallESB)) then
                f8(k+Nxtotal,l)=f8w(which_corner(ii,jj+1,kk),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(ii,jj,kk+1)==wallSB) .OR.(image(ii,jj,kk+1)==wallEB) &
            .OR.(image(ii,jj,kk+1)==wallENB).OR.(image(ii,jj,kk+1)==wallWSB)) then
                f8(k+Nxytotal,l)=f8w(which_corner(ii,jj,kk+1),l)
            else if (image(ii,jj,kk+1)==wallESB) then
                f8(k+Nxytotal,l)=f2(k+Nxytotal,l)
            End if

            fEq=w(l)*(Rho(k)+1.d0*(cx(l)*Ux(k)-cy(l)*Uy(k)-cz(l)*Uz(k)))
            f8(k,l)=(mu*(fEq-0.5d0*f8(k,l)) &
            &        + cx(l)*coef8(i,2)*f8(k-1,l) &
            &        + cx(l)*coef8(i,3)*f8(k-2,l) &
            &        - cy(l)*coef8(i,5)*f8(k+Nxtotal,l) &
            &        - cy(l)*coef8(i,6)*f8(k+2*Nxtotal,l) &
            &        - cz(l)*coef8(i,8)*f8(k+Nxytotal,l) &
            &        - cz(l)*coef8(i,9)*f8(k+2*Nxytotal,l) &
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

!$OMP DO 
    ! pack/unpack X dir buffer
    do k = 1,Nztotal
        do j = 1, Nytotal
            do i = 1, ghostLayers
                do l = 1, Nc8 ! dir (2,3,6,7) and (1,4,5,8)
                    ! pack order, 
                    !for dir: (2,3,6,7)
                    localidMsnd = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i+ghostLayers
                    localidMrcv = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i
                    !for dir: (1,4,5,8)
                    localidPsnd = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i+Nxsub
                    localidPrcv = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i+ghostLayers+Nxsub


                    packid = 0*Nztotal*Nytotal*ghostLayers*Nc8 & !dir 2, 1
                           + (k-1)*Nytotal*ghostLayers*Nc8 &
                           + (j-1)*ghostLayers*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_west_snd(packid) = f2(localidMsnd, l)
                    f2(localidPrcv,l) = f_east_rcv(packid)
                    f_east_snd(packid) = f1(localidPsnd, l)
                    f1(localidMrcv,l) = f_west_rcv(packid)

                    packid = 1*Nztotal*Nytotal*ghostLayers*Nc8 & !dir 3, 4
                           + (k-1)*Nytotal*ghostLayers*Nc8 &
                           + (j-1)*ghostLayers*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_west_snd(packid) = f3(localidMsnd, l)
                    f3(localidPrcv,l) = f_east_rcv(packid)
                    f_east_snd(packid) = f4(localidPsnd, l)
                    f4(localidMrcv,l) = f_west_rcv(packid)

                    packid = 2*Nztotal*Nytotal*ghostLayers*Nc8 & !dir 6, 5
                           + (k-1)*Nytotal*ghostLayers*Nc8 &
                           + (j-1)*ghostLayers*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_west_snd(packid) = f6(localidMsnd, l)
                    f6(localidPrcv,l) = f_east_rcv(packid)
                    f_east_snd(packid) = f5(localidPsnd, l)
                    f5(localidMrcv,l) = f_west_rcv(packid)

                    packid = 3*Nztotal*Nytotal*ghostLayers*Nc8 & !dir 7, 8
                           + (k-1)*Nytotal*ghostLayers*Nc8 &
                           + (j-1)*ghostLayers*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_west_snd(packid) = f7(localidMsnd, l)
                    f7(localidPrcv,l) = f_east_rcv(packid)
                    f_east_snd(packid) = f8(localidPsnd, l)
                    f8(localidMrcv,l) = f_west_rcv(packid)
                enddo !l
            enddo !i
        enddo !j
    enddo !k
!$OMP END DO

!$OMP DO 
    ! pack/unpack Y dir buffer
    do k = 1, Nztotal
        do j = 1, ghostLayers
            do i = 1, Nxtotal
                do l = 1, Nc8 ! dir (3,4,7,8) and (2,1,6,5)
                    ! pack order, 
                    !for dir: (3,4,7,8)
                    localidMsnd = (k-1)*Nxtotal*Nytotal + (j+ghostLayers-1)*Nxtotal + i
                    localidMrcv = (k-1)*Nxtotal*Nytotal + (j-1)*Nxtotal + i
                    !for dir: (2,1,6,5)
                    localidPsnd = (k-1)*Nxtotal*Nytotal + (j+Nysub-1)*Nxtotal + i
                    localidPrcv = (k-1)*Nxtotal*Nytotal + (j+Nysub+ghostLayers-1)*Nxtotal + i

                    packid = 0*Nztotal*ghostLayers*Nxtotal*Nc8 & !dir 3, 2
                           + (k-1)*ghostLayers*Nxtotal*Nc8 &
                           + (j-1)*Nxtotal*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_suth_snd(packid) = f3(localidMsnd, l)
                    f3(localidPrcv,l) = f_noth_rcv(packid)
                    f_noth_snd(packid) = f2(localidPsnd, l)
                    f2(localidMrcv,l) = f_suth_rcv(packid)

                    packid = 1*Nztotal*ghostLayers*Nxtotal*Nc8 & !dir 4, 1
                           + (k-1)*ghostLayers*Nxtotal*Nc8 &
                           + (j-1)*Nxtotal*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_suth_snd(packid) = f4(localidMsnd, l)
                    f4(localidPrcv,l) = f_noth_rcv(packid)
                    f_noth_snd(packid) = f1(localidPsnd, l)
                    f1(localidMrcv,l) = f_suth_rcv(packid)

                    packid = 2*Nztotal*ghostLayers*Nxtotal*Nc8 & !dir 7, 6
                           + (k-1)*ghostLayers*Nxtotal*Nc8 &
                           + (j-1)*Nxtotal*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_suth_snd(packid) = f7(localidMsnd, l)
                    f7(localidPrcv,l) = f_noth_rcv(packid)
                    f_noth_snd(packid) = f6(localidPsnd, l)
                    f6(localidMrcv,l) = f_suth_rcv(packid)


                    packid = 3*Nztotal*ghostLayers*Nxtotal*Nc8 & !dir 8, 5
                           + (k-1)*ghostLayers*Nxtotal*Nc8 &
                           + (j-1)*Nxtotal*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_suth_snd(packid) = f8(localidMsnd, l)
                    f8(localidPrcv,l) = f_noth_rcv(packid)
                    f_noth_snd(packid) = f5(localidPsnd, l)
                    f5(localidMrcv,l) = f_suth_rcv(packid)
                enddo !l
            enddo !i
        enddo !j
    enddo !k
!$OMP END DO


!$OMP DO 
    ! pack/unpack Z dir buffer
    do k = 1, ghostLayers
        do j = 1, Nytotal
            do i = 1, Nxtotal
                do l = 1, Nc8 ! dir (5,6,7,8) and (1,2,3,4)
                    ! pack order, 
                    !for dir: (5,6,7,8)
                    localidMsnd = (k+ghostLayers-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                    localidMrcv = (k-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                    !for dir: (1,2,3,4)
                    localidPsnd = (k+Nzsub-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i
                    localidPrcv = (k+Nzsub+ghostLayers-1)*Nytotal*Nxtotal + (j-1)*Nxtotal + i

                    packid = 0*ghostLayers*Nytotal*Nxtotal*Nc8 & !dir 5, 1
                           + (k-1)*Nytotal*Nxtotal*Nc8 &
                           + (j-1)*Nxtotal*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_back_snd(packid) = f5(localidMsnd, l)
                    f5(localidPrcv,l) = f_frnt_rcv(packid)
                    f_frnt_snd(packid) = f1(localidPsnd, l)
                    f1(localidMrcv,l) = f_back_rcv(packid)

                    packid = 1*ghostLayers*Nytotal*Nxtotal*Nc8 & !dir 6, 2
                           + (k-1)*Nytotal*Nxtotal*Nc8 &
                           + (j-1)*Nxtotal*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_back_snd(packid) = f6(localidMsnd, l)
                    f6(localidPrcv,l) = f_frnt_rcv(packid)
                    f_frnt_snd(packid) = f2(localidPsnd, l)
                    f2(localidMrcv,l) = f_back_rcv(packid) 

                    packid = 2*ghostLayers*Nytotal*Nxtotal*Nc8 & !dir 7, 3
                           + (k-1)*Nytotal*Nxtotal*Nc8 &
                           + (j-1)*Nxtotal*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_back_snd(packid) = f7(localidMsnd, l)
                    f7(localidPrcv,l) = f_frnt_rcv(packid)
                    f_frnt_snd(packid) = f3(localidPsnd, l)
                    f3(localidMrcv,l) = f_back_rcv(packid) 

                    packid = 3*ghostLayers*Nytotal*Nxtotal*Nc8 & !dir 8, 4
                           + (k-1)*Nytotal*Nxtotal*Nc8 &
                           + (j-1)*Nxtotal*Nc8 &
                           + (i-1)*Nc8 &
                           + l
                    f_back_snd(packid) = f8(localidMsnd, l)
                    f8(localidPrcv,l) = f_frnt_rcv(packid)
                    f_frnt_snd(packid) = f4(localidPsnd, l)
                    f4(localidMrcv,l) = f_back_rcv(packid)                   
                enddo !l
            enddo !i
        enddo !j
    enddo !k
!$OMP END DO


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
                    f2(k,l)=2.d0*f2(k+1,l)-f2(k+2,l)
                    f3(k,l)=2.d0*f3(k+1,l)-f3(k+2,l)
                    f6(k,l)=2.d0*f6(k+1,l)-f6(k+2,l)
                    f7(k,l)=2.d0*f7(k+1,l)-f7(k+2,l)
                    RhoWall=RhoWall+cx(l)*(f2(k,l)+f3(k,l)+f6(k,l)+f7(k,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                Enddo
            CASE (wallW)
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                    f4(k,l)=2.d0*f4(k-1,l)-f4(k-2,l)
                    f5(k,l)=2.d0*f5(k-1,l)-f5(k-2,l)
                    f8(k,l)=2.d0*f8(k-1,l)-f8(k-2,l)
                    RhoWall=RhoWall+cx(l)*(f1(k,l)+f4(k,l)+f5(k,l)+f8(k,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                Enddo
            CASE (wallN)
                Do l=1,Nc8
                    f3(k,l)=2.d0*f3(k+Nxtotal,l)-f3(k+2*Nxtotal,l)
                    f4(k,l)=2.d0*f4(k+Nxtotal,l)-f4(k+2*Nxtotal,l)
                    f7(k,l)=2.d0*f7(k+Nxtotal,l)-f7(k+2*Nxtotal,l)
                    f8(k,l)=2.d0*f8(k+Nxtotal,l)-f8(k+2*Nxtotal,l)
                    RhoWall=RhoWall+cy(l)*(f3(k,l)+f4(k,l)+f7(k,l)+f8(k,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                Enddo
            CASE (wallS)
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                    f2(k,l)=2.d0*f2(k-Nxtotal,l)-f2(k-2*Nxtotal,l)
                    f5(k,l)=2.d0*f5(k-Nxtotal,l)-f5(k-2*Nxtotal,l)
                    f6(k,l)=2.d0*f6(k-Nxtotal,l)-f6(k-2*Nxtotal,l)
                    RhoWall=RhoWall+cy(l)*(f1(k,l)+f2(k,l)+f5(k,l)+f6(k,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                Enddo
            CASE (wallF)
                Do l=1,Nc8
                    f5(k,l)=2.d0*f5(k+Nxytotal,l)-f5(k+2*Nxytotal,l)
                    f6(k,l)=2.d0*f6(k+Nxytotal,l)-f6(k+2*Nxytotal,l)
                    f7(k,l)=2.d0*f7(k+Nxytotal,l)-f7(k+2*Nxytotal,l)
                    f8(k,l)=2.d0*f8(k+Nxytotal,l)-f8(k+2*Nxytotal,l)
                    RhoWall=RhoWall+cz(l)*(f5(k,l)+f6(k,l)+f7(k,l)+f8(k,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                Enddo
            CASE (wallB)
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-Nxytotal,l)-f1(k-2*Nxytotal,l)
                    f2(k,l)=2.d0*f2(k-Nxytotal,l)-f2(k-2*Nxytotal,l)
                    f3(k,l)=2.d0*f3(k-Nxytotal,l)-f3(k-2*Nxytotal,l)
                    f4(k,l)=2.d0*f4(k-Nxytotal,l)-f4(k-2*Nxytotal,l)
                    RhoWall=RhoWall+cz(l)*(f1(k,l)+f2(k,l)+f3(k,l)+f4(k,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                Do l=1,Nc8
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                Enddo
!=======================================================================
!     Boundary condition on the 2-direction corner wall
!=======================================================================
!------------------------------------------------------------------------
!           No z-direction type wall
!------------------------------------------------------------------------
            CASE (wallEN)
                Do l=1,Nc8
                    f2(k,l)=2.d0*f2(k+1,l)-f2(k+2,l)
                    f3(k,l)=2.d0*f3(k+1,l)-f3(k+2,l)
                    f6(k,l)=2.d0*f6(k+1,l)-f6(k+2,l)
                    f7(k,l)=2.d0*f7(k+1,l)-f7(k+2,l)

                    f3w(j,l)=2.d0*f3(k+Nxtotal,l)-f3(k+2*Nxtotal,l)
                    f4w(j,l)=2.d0*f4(k+Nxtotal,l)-f4(k+2*Nxtotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxtotal,l)-f7(k+2*Nxtotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxtotal,l)-f8(k+2*Nxtotal,l)

                    RhoWall=RhoWall+cx(l)*(f2(k,l)+f3(k,l)+f6(k,l)+f7(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f3w(j,l)+f4w(j,l)+f7w(j,l)+f8w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write y
                    f2(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                    f6(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                    ! Store y
                    f1w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                    f5w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                Enddo

            CASE (wallWN)
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                    f4(k,l)=2.d0*f4(k-1,l)-f4(k-2,l)
                    f5(k,l)=2.d0*f5(k-1,l)-f5(k-2,l)
                    f8(k,l)=2.d0*f8(k-1,l)-f8(k-2,l)

                    f3w(j,l)=2.d0*f3(k+Nxtotal,l)-f3(k+2*Nxtotal,l)
                    f4w(j,l)=2.d0*f4(k+Nxtotal,l)-f4(k+2*Nxtotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxtotal,l)-f7(k+2*Nxtotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxtotal,l)-f8(k+2*Nxtotal,l)

                    RhoWall=RhoWall+cx(l)*(f1(k,l)+f4(k,l)+f5(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f3w(j,l)+f4w(j,l)+f7w(j,l)+f8w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8

                    ! Write x
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    ! Write y
                    f1(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                    f5(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                    ! Store y
                    f2w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                    f6w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                Enddo

            CASE (wallES)
                Do l=1,Nc8
                    f2(k,l)=2.d0*f2(k+1,l)-f2(k+2,l)
                    f3(k,l)=2.d0*f3(k+1,l)-f3(k+2,l)
                    f6(k,l)=2.d0*f6(k+1,l)-f6(k+2,l)
                    f7(k,l)=2.d0*f7(k+1,l)-f7(k+2,l)

                    f1w(j,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxtotal,l)-f2(k-2*Nxtotal,l)
                    f5w(j,l)=2.d0*f5(k-Nxtotal,l)-f5(k-2*Nxtotal,l)
                    f6w(j,l)=2.d0*f6(k-Nxtotal,l)-f6(k-2*Nxtotal,l)

                    RhoWall=RhoWall+cx(l)*(f2(k,l)+f3(k,l)+f6(k,l)+f7(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f1w(j,l)+f2w(j,l)+f5w(j,l)+f6w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write y
                    f3(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                    f7(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                    ! Store y
                    f4w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f8w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                Enddo

            CASE (wallWS)
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                    f4(k,l)=2.d0*f4(k-1,l)-f4(k-2,l)
                    f5(k,l)=2.d0*f5(k-1,l)-f5(k-2,l)
                    f8(k,l)=2.d0*f8(k-1,l)-f8(k-2,l)

                    f1w(j,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxtotal,l)-f2(k-2*Nxtotal,l)
                    f5w(j,l)=2.d0*f5(k-Nxtotal,l)-f5(k-2*Nxtotal,l)
                    f6w(j,l)=2.d0*f6(k-Nxtotal,l)-f6(k-2*Nxtotal,l)

                    RhoWall=RhoWall+cx(l)*(f1(k,l)+f4(k,l)+f5(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f1w(j,l)+f2w(j,l)+f5w(j,l)+f6w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    ! Write y
                    f4(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f8(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                    ! Store y
                    f3w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                    f7w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                Enddo

!------------------------------------------------------------------------
!           No y-direction type wall
!------------------------------------------------------------------------
            CASE (wallEF)
                Do l=1,Nc8
                    f2(k,l)=2.d0*f2(k+1,l)-f2(k+2,l)
                    f3(k,l)=2.d0*f3(k+1,l)-f3(k+2,l)
                    f6(k,l)=2.d0*f6(k+1,l)-f6(k+2,l)
                    f7(k,l)=2.d0*f7(k+1,l)-f7(k+2,l)

                    f5w(j,l)=2.d0*f5(k+Nxytotal,l)-f5(k+2*Nxytotal,l)
                    f6w(j,l)=2.d0*f6(k+Nxytotal,l)-f6(k+2*Nxytotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxytotal,l)-f7(k+2*Nxytotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxytotal,l)-f8(k+2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f2(k,l)+f3(k,l)+f6(k,l)+f7(k,l))
                    RhoWall2=RhoWall2+cz(l)*(f5w(j,l)+f6w(j,l)+f7w(j,l)+f8w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write z
                    f2(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                    f3(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                    ! Store z
                    f1w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                    f4w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                Enddo

            CASE (wallWF)
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                    f4(k,l)=2.d0*f4(k-1,l)-f4(k-2,l)
                    f5(k,l)=2.d0*f5(k-1,l)-f5(k-2,l)
                    f8(k,l)=2.d0*f8(k-1,l)-f8(k-2,l)

                    f5w(j,l)=2.d0*f5(k+Nxytotal,l)-f5(k+2*Nxytotal,l)
                    f6w(j,l)=2.d0*f6(k+Nxytotal,l)-f6(k+2*Nxytotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxytotal,l)-f7(k+2*Nxytotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxytotal,l)-f8(k+2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f1(k,l)+f4(k,l)+f5(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cz(l)*(f5w(j,l)+f6w(j,l)+f7w(j,l)+f8w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    ! Write z
                    f1(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                    f4(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                    ! Store z
                    f2w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                    f3w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                Enddo

            CASE (wallEB)
                Do l=1,Nc8
                    f2(k,l)=2.d0*f2(k+1,l)-f2(k+2,l)
                    f3(k,l)=2.d0*f3(k+1,l)-f3(k+2,l)
                    f6(k,l)=2.d0*f6(k+1,l)-f6(k+2,l)
                    f7(k,l)=2.d0*f7(k+1,l)-f7(k+2,l)

                    f1w(j,l)=2.d0*f1(k-Nxytotal,l)-f1(k-2*Nxytotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxytotal,l)-f2(k-2*Nxytotal,l)
                    f3w(j,l)=2.d0*f3(k-Nxytotal,l)-f3(k-2*Nxytotal,l)
                    f4w(j,l)=2.d0*f4(k-Nxytotal,l)-f4(k-2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f2(k,l)+f3(k,l)+f6(k,l)+f7(k,l))
                    RhoWall2=RhoWall2+cz(l)*(f1w(j,l)+f2w(j,l)+f3w(j,l)+f4w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write z
                    f6(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                    f7(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                    ! Store z
                    f5w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f8w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                Enddo

            CASE (wallWB)
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                    f4(k,l)=2.d0*f4(k-1,l)-f4(k-2,l)
                    f5(k,l)=2.d0*f5(k-1,l)-f5(k-2,l)
                    f8(k,l)=2.d0*f8(k-1,l)-f8(k-2,l)

                    f1w(j,l)=2.d0*f1(k-Nxytotal,l)-f1(k-2*Nxytotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxytotal,l)-f2(k-2*Nxytotal,l)
                    f3w(j,l)=2.d0*f3(k-Nxytotal,l)-f3(k-2*Nxytotal,l)
                    f4w(j,l)=2.d0*f4(k-Nxytotal,l)-f4(k-2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f1(k,l)+f4(k,l)+f5(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cz(l)*(f1w(j,l)+f2w(j,l)+f3w(j,l)+f4w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    ! Write z
                    f5(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f8(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                    ! Store z
                     f6w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                    f7w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                Enddo
!------------------------------------------------------------------------
!           No x-direction type wall
!------------------------------------------------------------------------
            CASE (wallNF)
                Do l=1,Nc8
                    f3(k,l)=2.d0*f3(k+Nxtotal,l)-f3(k+2*Nxtotal,l)
                    f4(k,l)=2.d0*f4(k+Nxtotal,l)-f4(k+2*Nxtotal,l)
                    f7(k,l)=2.d0*f7(k+Nxtotal,l)-f7(k+2*Nxtotal,l)
                    f8(k,l)=2.d0*f8(k+Nxtotal,l)-f8(k+2*Nxtotal,l)

                    f5w(j,l)=2.d0*f5(k+Nxytotal,l)-f5(k+2*Nxytotal,l)
                    f6w(j,l)=2.d0*f6(k+Nxytotal,l)-f6(k+2*Nxytotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxytotal,l)-f7(k+2*Nxytotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxytotal,l)-f8(k+2*Nxytotal,l)

                    RhoWall=RhoWall+cy(l)*(f3(k,l)+f4(k,l)+f7(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cz(l)*(f5w(j,l)+f6w(j,l)+f7w(j,l)+f8w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write y
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write z
                    f3(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                    f4(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                    ! Store z
                     f1w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                    f2w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                Enddo

            CASE (wallSF)
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                    f2(k,l)=2.d0*f2(k-Nxtotal,l)-f2(k-2*Nxtotal,l)
                    f5(k,l)=2.d0*f5(k-Nxtotal,l)-f5(k-2*Nxtotal,l)
                    f6(k,l)=2.d0*f6(k-Nxtotal,l)-f6(k-2*Nxtotal,l)

                    f5w(j,l)=2.d0*f5(k+Nxytotal,l)-f5(k+2*Nxytotal,l)
                    f6w(j,l)=2.d0*f6(k+Nxytotal,l)-f6(k+2*Nxytotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxytotal,l)-f7(k+2*Nxytotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxytotal,l)-f8(k+2*Nxytotal,l)

                    RhoWall=RhoWall+cy(l)*(f1(k,l)+f2(k,l)+f5(k,l)+f6(k,l))
                    RhoWall2=RhoWall2+cz(l)*(f5w(j,l)+f6w(j,l)+f7w(j,l)+f8w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write y
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    ! Write z
                    f1(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                    f2(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                    ! Store z
                    f3w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                    f4w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                Enddo

            CASE (wallNB)
                Do l=1,Nc8
                    f3(k,l)=2.d0*f3(k+Nxtotal,l)-f3(k+2*Nxtotal,l)
                    f4(k,l)=2.d0*f4(k+Nxtotal,l)-f4(k+2*Nxtotal,l)
                    f7(k,l)=2.d0*f7(k+Nxtotal,l)-f7(k+2*Nxtotal,l)
                    f8(k,l)=2.d0*f8(k+Nxtotal,l)-f8(k+2*Nxtotal,l)

                    f1w(j,l)=2.d0*f1(k-Nxytotal,l)-f1(k-2*Nxytotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxytotal,l)-f2(k-2*Nxytotal,l)
                    f3w(j,l)=2.d0*f3(k-Nxytotal,l)-f3(k-2*Nxytotal,l)
                    f4w(j,l)=2.d0*f4(k-Nxytotal,l)-f4(k-2*Nxytotal,l)

                    RhoWall=RhoWall+cy(l)*(f3(k,l)+f4(k,l)+f7(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cz(l)*(f1w(j,l)+f2w(j,l)+f3w(j,l)+f4w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write y
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write z
                    f7(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                    f8(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                    ! Store z
                    f5w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f6w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                Enddo

            CASE (wallSB)
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                    f2(k,l)=2.d0*f2(k-Nxtotal,l)-f2(k-2*Nxtotal,l)
                    f5(k,l)=2.d0*f5(k-Nxtotal,l)-f5(k-2*Nxtotal,l)
                    f6(k,l)=2.d0*f6(k-Nxtotal,l)-f6(k-2*Nxtotal,l)

                    f1w(j,l)=2.d0*f1(k-Nxytotal,l)-f1(k-2*Nxytotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxytotal,l)-f2(k-2*Nxytotal,l)
                    f3w(j,l)=2.d0*f3(k-Nxytotal,l)-f3(k-2*Nxytotal,l)
                    f4w(j,l)=2.d0*f4(k-Nxytotal,l)-f4(k-2*Nxytotal,l)

                    RhoWall=RhoWall+cy(l)*(f1(k,l)+f2(k,l)+f5(k,l)+f6(k,l))
                    RhoWall2=RhoWall2+cz(l)*(f1w(j,l)+f2w(j,l)+f3w(j,l)+f4w(j,l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                Do l=1,Nc8
                    ! Write y
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    ! Write z
                    f5(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f6(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                    ! Store z
                    f7w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                    f8w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                Enddo
!=======================================================================
!     Boundary condition on the 3-direction corner wall
!=======================================================================
            CASE (wallENF) !direction1
                Do l=1,Nc8
                    f2(k,l)=2.d0*f2(k+1,l)-f2(k+2,l)
                    f3(k,l)=2.d0*f3(k+1,l)-f3(k+2,l)
                    f6(k,l)=2.d0*f6(k+1,l)-f6(k+2,l)
                    f7(k,l)=2.d0*f7(k+1,l)-f7(k+2,l)

                    f3w(j,l)=2.d0*f3(k+Nxtotal,l)-f3(k+2*Nxtotal,l)
                    f4w(j,l)=2.d0*f4(k+Nxtotal,l)-f4(k+2*Nxtotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxtotal,l)-f7(k+2*Nxtotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxtotal,l)-f8(k+2*Nxtotal,l)

                    f5wZ(l)=2.d0*f5(k+Nxytotal,l)-f5(k+2*Nxytotal,l)
                    f6wZ(l)=2.d0*f6(k+Nxytotal,l)-f6(k+2*Nxytotal,l)
                    f7wZ(l)=2.d0*f7(k+Nxytotal,l)-f7(k+2*Nxytotal,l)
                    f8wZ(l)=2.d0*f8(k+Nxytotal,l)-f8(k+2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f2(k,l)+f3(k,l)+f6(k,l)+f7(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f3w(j,l)+f4w(j,l)+f7w(j,l)+f8w(j,l))
                    RhoWall3=RhoWall3+cz(l)*(f5wZ(l)+f6wZ(l)+f7wZ(l)+f8wZ(l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write y
                    f2(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                    f6(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                    ! Write z
                    f3(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f7wZ(l)

                    ! Store y
                    f1w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                    f5w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                    ! Store z
                    f2w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f6wZ(l)
                    f4w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f8wZ(l)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f7(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f5wZ(l)
                Enddo
            CASE (wallWNF) !direction2
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                    f4(k,l)=2.d0*f4(k-1,l)-f4(k-2,l)
                    f5(k,l)=2.d0*f5(k-1,l)-f5(k-2,l)
                    f8(k,l)=2.d0*f8(k-1,l)-f8(k-2,l)

                    f3w(j,l)=2.d0*f3(k+Nxtotal,l)-f3(k+2*Nxtotal,l)
                    f4w(j,l)=2.d0*f4(k+Nxtotal,l)-f4(k+2*Nxtotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxtotal,l)-f7(k+2*Nxtotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxtotal,l)-f8(k+2*Nxtotal,l)

                    f5wZ(l)=2.d0*f5(k+Nxytotal,l)-f5(k+2*Nxytotal,l)
                    f6wZ(l)=2.d0*f6(k+Nxytotal,l)-f6(k+2*Nxytotal,l)
                    f7wZ(l)=2.d0*f7(k+Nxytotal,l)-f7(k+2*Nxytotal,l)
                    f8wZ(l)=2.d0*f8(k+Nxytotal,l)-f8(k+2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f1(k,l)+f4(k,l)+f5(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f3w(j,l)+f4w(j,l)+f7w(j,l)+f8w(j,l))
                    RhoWall3=RhoWall3+cz(l)*(f5wZ(l)+f6wZ(l)+f7wZ(l)+f8wZ(l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    ! Write y
                    f1(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                    f5(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                    ! Write z
                    f4(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f8wZ(l)

                    ! Store y
                    f2w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                    f6w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                    ! Store z
                    f1w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f5wZ(l)
                    f3w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f7wZ(l)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f8(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f6wZ(l)
                Enddo
            CASE (wallWSF) !direction3
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                    f4(k,l)=2.d0*f4(k-1,l)-f4(k-2,l)
                    f5(k,l)=2.d0*f5(k-1,l)-f5(k-2,l)
                    f8(k,l)=2.d0*f8(k-1,l)-f8(k-2,l)

                    f1w(j,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxtotal,l)-f2(k-2*Nxtotal,l)
                    f5w(j,l)=2.d0*f5(k-Nxtotal,l)-f5(k-2*Nxtotal,l)
                    f6w(j,l)=2.d0*f6(k-Nxtotal,l)-f6(k-2*Nxtotal,l)

                    f5wZ(l)=2.d0*f5(k+Nxytotal,l)-f5(k+2*Nxytotal,l)
                    f6wZ(l)=2.d0*f6(k+Nxytotal,l)-f6(k+2*Nxytotal,l)
                    f7wZ(l)=2.d0*f7(k+Nxytotal,l)-f7(k+2*Nxytotal,l)
                    f8wZ(l)=2.d0*f8(k+Nxytotal,l)-f8(k+2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f1(k,l)+f4(k,l)+f5(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f1w(j,l)+f2w(j,l)+f5w(j,l)+f6w(j,l))
                    RhoWall3=RhoWall3+cz(l)*(f5wZ(l)+f6wZ(l)+f7wZ(l)+f8wZ(l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    ! Write y
                    f4(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f8(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                    ! Write z
                    f1(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f5wZ(l)

                    ! Store y
                    f3w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                    f7w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                    ! Store z
                    f2w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f6wZ(l)
                    f4w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f8wZ(l)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f5(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f7wZ(l)
                Enddo
            CASE (wallESF) !direction4
                Do l=1,Nc8
                    f2(k,l)=2.d0*f2(k+1,l)-f2(k+2,l)
                    f3(k,l)=2.d0*f3(k+1,l)-f3(k+2,l)
                    f6(k,l)=2.d0*f6(k+1,l)-f6(k+2,l)
                    f7(k,l)=2.d0*f7(k+1,l)-f7(k+2,l)

                    f1w(j,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxtotal,l)-f2(k-2*Nxtotal,l)
                    f5w(j,l)=2.d0*f5(k-Nxtotal,l)-f5(k-2*Nxtotal,l)
                    f6w(j,l)=2.d0*f6(k-Nxtotal,l)-f6(k-2*Nxtotal,l)

                    f5wZ(l)=2.d0*f5(k+Nxytotal,l)-f5(k+2*Nxytotal,l)
                    f6wZ(l)=2.d0*f6(k+Nxytotal,l)-f6(k+2*Nxytotal,l)
                    f7wZ(l)=2.d0*f7(k+Nxytotal,l)-f7(k+2*Nxytotal,l)
                    f8wZ(l)=2.d0*f8(k+Nxytotal,l)-f8(k+2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f2(k,l)+f3(k,l)+f6(k,l)+f7(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f1w(j,l)+f2w(j,l)+f5w(j,l)+f6w(j,l))
                    RhoWall3=RhoWall3+cz(l)*(f5wZ(l)+f6wZ(l)+f7wZ(l)+f8wZ(l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write y
                    f3(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                    f7(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                    ! Write z
                    f2(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f6wZ(l)

                    ! Store y
                    f4w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f8w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                    ! Store z
                    f1w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f5wZ(l)
                    f3w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f7wZ(l)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f6(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f8wZ(l)
                Enddo
            CASE (wallENB) !direction5
                Do l=1,Nc8
                    f2(k,l)=2.d0*f2(k+1,l)-f2(k+2,l)
                    f3(k,l)=2.d0*f3(k+1,l)-f3(k+2,l)
                    f6(k,l)=2.d0*f6(k+1,l)-f6(k+2,l)
                    f7(k,l)=2.d0*f7(k+1,l)-f7(k+2,l)

                    f3w(j,l)=2.d0*f3(k+Nxtotal,l)-f3(k+2*Nxtotal,l)
                    f4w(j,l)=2.d0*f4(k+Nxtotal,l)-f4(k+2*Nxtotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxtotal,l)-f7(k+2*Nxtotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxtotal,l)-f8(k+2*Nxtotal,l)

                    f1wZ(l)=2.d0*f1(k-Nxytotal,l)-f1(k-2*Nxytotal,l)
                    f2wZ(l)=2.d0*f2(k-Nxytotal,l)-f2(k-2*Nxytotal,l)
                    f3wZ(l)=2.d0*f3(k-Nxytotal,l)-f3(k-2*Nxytotal,l)
                    f4wZ(l)=2.d0*f4(k-Nxytotal,l)-f4(k-2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f2(k,l)+f3(k,l)+f6(k,l)+f7(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f3w(j,l)+f4w(j,l)+f7w(j,l)+f8w(j,l))
                    RhoWall3=RhoWall3+cz(l)*(f1wZ(l)+f2wZ(l)+f3wZ(l)+f4wZ(l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write y
                    f2(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                    f6(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                    ! Write z
                    f7(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f3wZ(l)

                    ! Store y
                    f1w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                    f5w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                    ! Store z
                    f6w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f2wZ(l)
                    f8w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f4wZ(l)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f3(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f1wZ(l)
                Enddo
            CASE (wallWNB) !direction6
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                    f4(k,l)=2.d0*f4(k-1,l)-f4(k-2,l)
                    f5(k,l)=2.d0*f5(k-1,l)-f5(k-2,l)
                    f8(k,l)=2.d0*f8(k-1,l)-f8(k-2,l)

                    f3w(j,l)=2.d0*f3(k+Nxtotal,l)-f3(k+2*Nxtotal,l)
                    f4w(j,l)=2.d0*f4(k+Nxtotal,l)-f4(k+2*Nxtotal,l)
                    f7w(j,l)=2.d0*f7(k+Nxtotal,l)-f7(k+2*Nxtotal,l)
                    f8w(j,l)=2.d0*f8(k+Nxtotal,l)-f8(k+2*Nxtotal,l)

                    f1wZ(l)=2.d0*f1(k-Nxytotal,l)-f1(k-2*Nxytotal,l)
                    f2wZ(l)=2.d0*f2(k-Nxytotal,l)-f2(k-2*Nxytotal,l)
                    f3wZ(l)=2.d0*f3(k-Nxytotal,l)-f3(k-2*Nxytotal,l)
                    f4wZ(l)=2.d0*f4(k-Nxytotal,l)-f4(k-2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f1(k,l)+f4(k,l)+f5(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f3w(j,l)+f4w(j,l)+f7w(j,l)+f8w(j,l))
                    RhoWall3=RhoWall3+cz(l)*(f1wZ(l)+f2wZ(l)+f3wZ(l)+f4wZ(l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    ! Write y
                    f1(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f4w(j,l)
                    f5(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f8w(j,l)
                    ! Write z
                    f8(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f4wZ(l)

                    ! Store y
                    f2w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f3w(j,l)
                    f6w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f7w(j,l)
                    ! Store z
                    f5w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f1wZ(l)
                    f7w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f3wZ(l)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f4(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f2wZ(l)
                Enddo
            CASE (wallWSB) !direction7
                Do l=1,Nc8
                    f1(k,l)=2.d0*f1(k-1,l)-f1(k-2,l)
                    f4(k,l)=2.d0*f4(k-1,l)-f4(k-2,l)
                    f5(k,l)=2.d0*f5(k-1,l)-f5(k-2,l)
                    f8(k,l)=2.d0*f8(k-1,l)-f8(k-2,l)

                    f1w(j,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxtotal,l)-f2(k-2*Nxtotal,l)
                    f5w(j,l)=2.d0*f5(k-Nxtotal,l)-f5(k-2*Nxtotal,l)
                    f6w(j,l)=2.d0*f6(k-Nxtotal,l)-f6(k-2*Nxtotal,l)

                    f1wZ(l)=2.d0*f1(k-Nxytotal,l)-f1(k-2*Nxytotal,l)
                    f2wZ(l)=2.d0*f2(k-Nxytotal,l)-f2(k-2*Nxytotal,l)
                    f3wZ(l)=2.d0*f3(k-Nxytotal,l)-f3(k-2*Nxytotal,l)
                    f4wZ(l)=2.d0*f4(k-Nxytotal,l)-f4(k-2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f1(k,l)+f4(k,l)+f5(k,l)+f8(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f1w(j,l)+f2w(j,l)+f5w(j,l)+f6w(j,l))
                    RhoWall3=RhoWall3+cz(l)*(f1wZ(l)+f2wZ(l)+f3wZ(l)+f4wZ(l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f2(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f1(k,l)
                    f3(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f4(k,l)
                    f6(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f5(k,l)
                    f7(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f8(k,l)
                    ! Write y
                    f4(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f8(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                    ! Write z
                    f5(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f1wZ(l)

                    ! Store y
                    f3w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                    f7w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                    ! Store z
                    f6w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f2wZ(l)
                    f8w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f4wZ(l)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f1(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f3wZ(l)
                Enddo
            CASE (wallESB) !direction8
                Do l=1,Nc8
                    f2(k,l)=2.d0*f2(k+1,l)-f2(k+2,l)
                    f3(k,l)=2.d0*f3(k+1,l)-f3(k+2,l)
                    f6(k,l)=2.d0*f6(k+1,l)-f6(k+2,l)
                    f7(k,l)=2.d0*f7(k+1,l)-f7(k+2,l)

                    f1w(j,l)=2.d0*f1(k-Nxtotal,l)-f1(k-2*Nxtotal,l)
                    f2w(j,l)=2.d0*f2(k-Nxtotal,l)-f2(k-2*Nxtotal,l)
                    f5w(j,l)=2.d0*f5(k-Nxtotal,l)-f5(k-2*Nxtotal,l)
                    f6w(j,l)=2.d0*f6(k-Nxtotal,l)-f6(k-2*Nxtotal,l)

                    f1wZ(l)=2.d0*f1(k-Nxytotal,l)-f1(k-2*Nxytotal,l)
                    f2wZ(l)=2.d0*f2(k-Nxytotal,l)-f2(k-2*Nxytotal,l)
                    f3wZ(l)=2.d0*f3(k-Nxytotal,l)-f3(k-2*Nxytotal,l)
                    f4wZ(l)=2.d0*f4(k-Nxytotal,l)-f4(k-2*Nxytotal,l)

                    RhoWall=RhoWall+cx(l)*(f2(k,l)+f3(k,l)+f6(k,l)+f7(k,l))
                    RhoWall2=RhoWall2+cy(l)*(f1w(j,l)+f2w(j,l)+f5w(j,l)+f6w(j,l))
                    RhoWall3=RhoWall3+cz(l)*(f1wZ(l)+f2wZ(l)+f3wZ(l)+f4wZ(l))
                Enddo
                RhoWall=RhoWall/DiffFlux
                RhoWall2=RhoWall2/DiffFlux
                RhoWall3=RhoWall3/DiffFlux
                Do l=1,Nc8
                    ! Write x
                    f1(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f2(k,l)
                    f4(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f3(k,l)
                    f5(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f6(k,l)
                    f8(k,l)=accom*w(l)*RhoWall &
                    &       + (1.d0-accom)*f7(k,l)
                    ! Write y
                    f3(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f2w(j,l)
                    f7(k,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f6w(j,l)
                    ! Write z
                    f6(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f2wZ(l)

                    ! Store y
                    f4w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f1w(j,l)
                    f8w(j,l)=accom*w(l)*RhoWall2 &
                    &       + (1.d0-accom)*f5w(j,l)
                    ! Store z
                    f5w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f1wZ(l)
                    f7w(j,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f3wZ(l)
                    ! Borrow the un-used variable to store z for the opposite velocity
                    f2(k,l)=accom*w(l)*RhoWall3 &
                    &       + (1.d0-accom)*f4wZ(l)
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
                f1(l,:)=f1(l-1,:) + w(:)*PressDrop
                f4(l,:)=f4(l-1,:) + w(:)*PressDrop
                f5(l,:)=f5(l-1,:) + w(:)*PressDrop
                f8(l,:)=f8(l-1,:) + w(:)*PressDrop
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
                f2(l,:)=f2(l+1,:) - w(:)*PressDrop
                f3(l,:)=f3(l+1,:) - w(:)*PressDrop
                f6(l,:)=f6(l+1,:) - w(:)*PressDrop
                f7(l,:)=f7(l+1,:) - w(:)*PressDrop
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
                f2(l,:)=f3(l,:)
                f1(l,:)=f4(l,:)
                f6(l,:)=f7(l,:)
                f5(l,:)=f8(l,:)
            enddo
        End do
!$OMP END DO 
    endif
    if(yu==ymax) then ! north sym BC
!$OMP DO
        Do k=zl,zu
            Do i=xl,xu
                l = (k-zlg)*Nxytotal + (ghostLayers+Nysub-1)*Nxtotal + i-xlg+1
                f3(l,:)=f2(l,:)
                f4(l,:)=f1(l,:)
                f7(l,:)=f6(l,:)
                f8(l,:)=f5(l,:)
            End do
        End do
!$OMP END DO
    endif
    if(zl==zmin) then ! back sym BC
!$OMP DO
        Do j=yl,yu
            Do i=xl,xu
                l = ghostLayers*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                f1(l,:)=f5(l,:)
                f2(l,:)=f6(l,:)
                f3(l,:)=f7(l,:)
                f4(l,:)=f8(l,:)
            enddo
        End do
!$OMP END DO 
    endif
    if(zu==zmax) then ! front sym BC
!$OMP DO
        Do j=yl,yu
            Do i=xl,xu
                l = (ghostLayers+Nzsub-1)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                f5(l,:)=f1(l,:)
                f6(l,:)=f2(l,:)
                f7(l,:)=f3(l,:)
                f8(l,:)=f4(l,:)
            End do
        End do
!$OMP END DO
    endif

!----------------------------------------------------
!> Update Macro
!----------------------------------------------------
!$OMP DO SCHEDULE(STATIC) 
    Do k=1,Ntotal
        Rho(k)=0.d0
        Ux(k)=0.d0
        Uy(k)=0.d0
        Uz(k)=0.d0
        Do l=1,Nc8
            Rho(k)=Rho(k)+f1(k,l)+f2(k,l)+f3(k,l)+f4(k,l)+f5(k,l)+f6(k,l)+f7(k,l)+f8(k,l)
            Ux(k)=Ux(k)+cx(l)*(f1(k,l)-f2(k,l)-f3(k,l)+f4(k,l)+f5(k,l)-f6(k,l)-f7(k,l)+f8(k,l))
            Uy(k)=Uy(k)+cy(l)*(f1(k,l)+f2(k,l)-f3(k,l)-f4(k,l)+f5(k,l)+f6(k,l)-f7(k,l)-f8(k,l))
            Uz(k)=Uz(k)+cz(l)*(f1(k,l)+f2(k,l)+f3(k,l)+f4(k,l)-f5(k,l)-f6(k,l)-f7(k,l)-f8(k,l))
        End do


    End do
!$OMP END DO 
!$OMP END PARALLEL  
    end subroutine iterate
   

    subroutine chkConverge
        implicit none
        include "mpif.h"


        ! local vars
        integer :: byl, byu, bzl, bzu, k, j, i, l, MPI_ERR
        double precision :: massInner, massNoth, massSuth, massFrnt, massBack
        double precision :: massNF, massNB, massSF, massSB
        double precision :: massLocal, mass2

        massInner = 0.d0
        massNoth = 0.d0
        massSuth = 0.d0
        massFrnt = 0.d0
        massBack = 0.d0
        massNF = 0.d0
        massNB = 0.d0
        massSF = 0.d0
        massSB = 0.d0
        
        byl = yl
        byu = yu
        bzl = zl
        bzu = zu
        if(yl == ymin) byl = yl + 1 !if most south block
        if(yu == ymax) byu = yu - 1 !if most north block
        if(zl == zmin) bzl = zl + 1 !if most back block
        if(zu == zmax) bzu = zu - 1 !if most frnt block

        !mass2=0.d0
        if(xl == xmin) then !only left most processors
            do k=bzl,bzu
                do j=byl,byu
                    l=(k-zlg)*Nxytotal + (j-ylg)*Nxtotal + column+ghostLayers
                    massInner=massInner+Ux(l)*ds*ds
                enddo
            enddo
            if (yl == ymin) then !only south most processors
                do k=bzl,bzu
                    l = (k-zlg)*Nxytotal + ghostLayers*Nxtotal + column+ghostLayers
                    massSuth = massSuth + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (yu == ymax) then !only north most processors
                do k=bzl,bzu
                    l = (k-zlg)*Nxytotal + (Nysub+ghostLayers-1)*Nxtotal + column+ghostLayers
                    massNoth = massNoth + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (zl == zmin) then !only back most processors
                do j=byl,byu
                    l = ghostLayers*Nxytotal + (j-ylg)*Nxtotal + column+ghostLayers
                    massBack = massBack + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (zu == zmax) then !only frnt most processors
                do j=byl,byu
                    l = (Nzsub+ghostLayers-1)*Nxytotal + (j-ylg)*Nxtotal + column+ghostLayers
                    massFrnt = massFrnt + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (yl==ymin.and.zl==zmin) then
                l = ghostLayers*Nxytotal + ghostLayers*Nxtotal + column+ghostLayers
                massSB = Ux(l)*ds*ds*0.25d0
            endif
            if (yl==ymin.and.zu==zmax) then
                l = (Nzsub+ghostLayers-1)*Nxytotal + ghostLayers*Nxtotal + column+ghostLayers
                massSF = Ux(l)*ds*ds*0.25d0
            endif
            if (yu==ymax.and.zl==zmin) then
                l = ghostLayers*Nxytotal + (Nysub+ghostLayers-1)*Nxtotal + column+ghostLayers
                massNB = Ux(l)*ds*ds*0.25d0
            endif
            if (yu==ymax.and.zu==zmax) then
                l = (Nzsub+ghostLayers-1)*Nxytotal + (Nysub+ghostLayers-1)*Nxtotal + column+ghostLayers
                massNF = Ux(l)*ds*ds*0.25d0
            endif


            ! debug
            massLocal = (massInner + massSuth + massNoth + massFrnt + massBack &
                + massSB + massSF + massNB + massNF) &
                *dsqrt(1.d0/2.d0)/PressDrop/(1.d0/Ref_L)**2
                !+ massSB + massSF + massNB + massNF)*dsqrt(1.d0/2.d0)*4.d0/1.d0

            ! reduction
            call MPI_ALLREDUCE(massLocal, mass2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               mpi_comm_inlet, MPI_ERR)

            !PRINT*, "mass = ", mass
            !DEBUG
            error=dabs(1.d0-mass2/mass)/(chkConvergeStep)

            mass=mass2
            if (proc == master) then           
                permeability=mass*Kn*sqrt(4.d0/PI)  
                write(*,"( 1I10, 3ES15.6)")  iStep,  mass,  permeability, error
                ! open(22,file='Results.dat', position="append")
                ! write(22,'(4ES15.6, 1I15)') Kn, mass, permeability, error, iStep
                ! close(22)
            endif
        endif ! xl==xmin
        !bcast error so every process in WORLD can stop
        CALL MPI_BCAST(error, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_VGRID, MPI_ERR)
    end subroutine chkConverge

    subroutine saveFlowField
        integer :: j, i, k, l
        character(13) fname
        character(9) ztitle
        write(fname, '(A, I0.3, A)') 'Field_', proc, '.dat'
        write(ztitle, '(A, I0.3,A)') '"proc', proc, '"'
        open(20,file=fname,STATUS="REPLACE")
        write(20,*) ' TITLE=" Field"'
        write(20,*) ' VARIABLES=x,y,z,flag,Rho,Ux,Uy,Uz'
        write(20,'(A,A,A,I,A,I,A,I,A)') ' ZONE T=', ztitle, ', I=', Nxsub,', J=', Nysub, ', K=', Nzsub, ' F=POINT'
        Do k=zl,zu
            Do j=yl,yu
                Do i=xl,xu
                    l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                    If (image(i,j,k)==fluid) then
                        write(20,'(4I10,4ES15.6)') i, j, k, 0, Rho(l)+1.d0, Ux(l), Uy(l), Uz(l)
                    else
                        write(20,'(4I10,4ES15.6)') i, j, k, 1, 0.d0, 0.d0, 0.d0, 0.d0
                    Endif   
                Enddo
            Enddo
        Enddo
        close(20)
    end subroutine saveFlowField


    SUBROUTINE saveFlowFieldVTK
        IMPLICIT NONE
        INTEGER :: i, j, k, l,MPI_ERR, IO_ERR
        character(13) fname

        write(fname, '(A, I0.3, A)') 'Field_', proc, '.vtk'
        OPEN(UNIT = 12, FILE = fname, STATUS = "REPLACE", POSITION = "APPEND", &
          IOSTAT = IO_ERR)
        IF ( IO_ERR == 0 ) THEN
            WRITE(12,'(A)')"# vtk DataFile Version 2.0"
            WRITE(12,'(A)')"DVM MPI"
            WRITE(12,'(A)')"ASCII"
            WRITE(12,*)
            WRITE(12,'(A)')"DATASET STRUCTURED_POINTS"
            WRITE(12,*)"DIMENSIONS",Nxsub,Nysub,Nzsub
            WRITE(12,*)"ORIGIN",xl,yl,zl
            WRITE(12,*)"SPACING",1,1,1
            WRITE(12,*)
            WRITE(12,*)"POINT_DATA",Nxsub*Nysub*Nzsub
            WRITE(12,*)
            WRITE(12,'(A)')"SCALARS Rho double"
            WRITE(12,'(A)')"LOOKUP_TABLE default"
            DO k = zl, zu
                DO j = yl, yu
                    DO i = xl, xu
                        l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                        If (image(i,j,k)==fluid) then
                            write(12,'(ES15.6)') Rho(l)+1.d0
                        else
                            write(12,'(ES15.6)') 0.d0
                        Endif   
                    END DO
                END DO
            END DO
            WRITE(12,*)
            WRITE(12,'(A)')"VECTORS Velocity double"
            DO k = zl, zu
               DO j = yl, yu
                  DO i = xl, xu
                     l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                     If (image(i,j,k)==fluid) then
                         write(12,'(3ES15.6)') Ux(l), Uy(l), Uz(l)
                     else
                         write(12,'(3ES15.6)') 0.d0, 0.d0, 0.d0
                     Endif  
                  END DO
               END DO
            END DO
            CLOSE(UNIT = 12)
        ELSE
            CALL memFree
            CALL MPI_FINALIZE(MPI_ERR)
            STOP "Error: Unable to open output vtk file."
        END IF

        RETURN
    END SUBROUTINE saveFlowFieldVTK

    SUBROUTINE saveFlowFieldVTI
        INTEGER :: i, j, k, l,MPI_ERR, IO_ERR
        character(13) fname
        integer :: exl, exu, eyl, eyu, ezl, ezu
        integer :: wexl, wexu, weyl, weyu, wezl, wezu
        wexl = xmin - 1
        wexu = xmax
        weyl = ymin - 1
        weyu = ymax
        wezl = zmin - 1
        wezu = zmax
        exl = xl - 1
        exu = xu
        eyl = yl - 1
        eyu = yu
        ezl = zl - 1
        ezu = zu

        write(fname, '(A, I0.3, A)') 'Field_', proc, '.vti'
        OPEN(UNIT = 13, FILE = fname, STATUS = "REPLACE", POSITION = "APPEND", &
          IOSTAT = IO_ERR)
        IF ( IO_ERR == 0 ) THEN
            WRITE(13,'(A)') '<?xml version="1.0"?>'
            WRITE(13,'(A)') '<VTKFile type="ImageData">'
            WRITE(13,'(A, 6I4, A)') '<ImageData WholeExtent="', exl, exu, & 
                eyl, eyu, ezl, ezu, ' " Origin="0 0 0" Spacing="1 1 1">'
            WRITE(13, '(A, 6I4, A)') '<Piece Extent="', exl, exu, eyl, eyu, &
                ezl, ezu, '">'
            WRITE(13, '(A)') '<CellData Scalars="flag Rho" Vectors ="U">'
            WRITE(13, '(A)') '<DataArray type="Int32" Name="flag" format="ascii">'
            DO k = zl, zu
                DO j = yl, yu
                    DO i = xl, xu
                        write(13,'(I)') image(i,j,k)
                    END DO
                END DO
            END DO
            WRITE(13, '(A)') '</DataArray>'
            WRITE(13, '(A)') '<DataArray type="Float32" Name="Rho" format="ascii">'
            DO k = zl, zu
                DO j = yl, yu
                    DO i = xl, xu
                        l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                        If (image(i,j,k)==fluid) then
                            write(13,'(ES15.6)') Rho(l)+1.d0
                        else
                            write(13,'(ES15.6)') 0.d0
                        Endif   
                    END DO
                END DO
            END DO
            WRITE(13, '(A)') '</DataArray>'
            WRITE(13, '(A)') '<DataArray type="Float32" Name="U" format="ascii" &
                & NumberOfComponents="3">'
            DO k = zl, zu
               DO j = yl, yu
                  DO i = xl, xu
                     l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                     If (image(i,j,k)==fluid) then
                         write(13,'(3ES15.6)') Ux(l), Uy(l), Uz(l)
                     else
                         write(13,'(3ES15.6)') 0.d0, 0.d0, 0.d0
                     Endif  
                  END DO
               END DO
            END DO
            WRITE(13, '(A)') '</DataArray>'
            WRITE(13, '(A)') '</CellData>'
            WRITE(13, '(A)') '</Piece>'
            WRITE(13, '(A)') '</ImageData>'
            WRITE(13, '(A)') '</VTKFile>'
            CLOSE(UNIT = 13)
        ELSE
            CALL memFree
            CALL MPI_FINALIZE(MPI_ERR)
            STOP "Error: Unable to open output vti file."
        END IF
    END SUBROUTINE saveFlowFieldVTI


    subroutine memFree
        deallocate (array3D, array3Dg)
        deallocate (image, which_corner, vecWall)
        deallocate (dir1, dir2, dir3, dir4, dir5, dir6, dir7, dir8)
        deallocate (coef1, coef2, coef3, coef4, coef5, coef6, coef7, coef8)
        deallocate (f1,f2,f3,f4,f5,f6,f7,f8)
        !deallocate (f1w,f2w,f3w,f4w,f5w,f6w,f7w,f8w)
        deallocate (Rho, Ux, Uy, Uz)
        deallocate (f_west_snd, f_west_rcv, f_east_snd, f_east_rcv)
        deallocate (f_noth_snd, f_noth_rcv, f_suth_snd, f_suth_rcv)
        deallocate (f_back_snd, f_back_rcv, f_frnt_snd, f_frnt_rcv)
    end subroutine memFree
end module solver
