!=======================================================================
!> @brief Physical space configurations
!=======================================================================
module physicalGrid
IMPLICIT NONE
SAVE
!--------------------------------------------------------------------
! Domain size
!domain defination, to be read from NML: velocityNml
!--------------------------------------------------------------------
! Total domain size (the raw image dimension with possible repeation)
integer :: Nx, Ny, Nz
! image file name
character(len=40):: imageFileName
! wall extrapolation order, 1, 2(default), 3(mixed)
integer :: wallExtOrder
! number of fluid layer to be added at inlet and outlet
integer :: fluidlayer
!1.d0 for reference length L=Lx=Ly or 2.d0 for L=Lx=2Ly
double precision:: Ref_L


! will be set by mpiParams
integer :: xl, xu, yl, yu, zl, zu, xlg, xug, ylg, yug, zlg, zug
! will be set in this module's init func.
integer :: xmin, xmax, ymin, ymax, zmin, zmax
integer :: block_repx, block_repy, block_repz
integer :: Nx_base, Ny_base, Nz_base

integer, parameter :: ghostLayers = 2
integer, parameter :: translate = 0
integer, parameter :: IMAGEFILE = 20

! need to be determined from mpi cood
integer :: Nxtotal, Nytotal, Nztotal, Nxytotal, Nxsub, Nysub, Nzsub, Ntotal
!uniform grid spacing in physical space, will be set in init func.
double precision :: ds
integer, parameter :: column= 2 ! layer to extract flow rate

! raw and extended flag array
integer, dimension (:,:,:), allocatable :: array3D, array3Dg ! globle image
integer, dimension (:,:,:), allocatable :: image, which_corner! local image 

! grid point flags and wall type id
integer, parameter:: fluid = 0, solid = 1, ghost=4
integer, parameter:: wallE = 20, wallW = 21, wallN = 22, wallS = 23, &
    wallF = 24, wallB = 25 !(1-side wall )
integer, parameter:: wallEN = 30, wallWN = 31, wallES = 32, wallWS = 33 !(2-side wall)
integer, parameter:: wallNF = 40, wallNB = 41, wallSF = 42, wallSB = 43 !(2-side wall)
integer, parameter:: wallEF = 50, wallEB = 51, wallWF = 52, wallWB = 53 !(2-side wall)
integer, parameter:: wallENF = 60, wallWNF = 61, wallWSF = 62, wallESF = 63 !(3-side wall)
integer, parameter:: wallENB = 70, wallWNB = 71, wallWSB = 72, wallESB = 73 !(3-side wall)
integer, parameter:: allWallLabels(26) = (/&
    wallE, wallW, wallN, wallS, wallF, wallB, &
    WallEN, wallWN, wallES, wallWS, &
    wallNF, wallNB, wallSF, wallSB, &
    wallEF, wallEB, wallWF, wallWB, &
    wallENF, wallWNF, wallWSF, wallESF, &
    wallENB, wallWNB, wallWSB, wallESB /)
integer, parameter:: onlyCornerLabels(20) = (/ &
    WallEN, wallWN, wallES, wallWS, &
    wallNF, wallNB, wallSF, wallSB, &
    wallEF, wallEB, wallWF, wallWB, &
    wallENF, wallWNF, wallWSF, wallESF, &
    wallENB, wallWNB, wallWSB, wallESB /)

integer :: Nstencil1, Nstencil2, Nstencil3, Nstencil4 !(group 1-4)
integer :: Nstencil5, Nstencil6, Nstencil7, Nstencil8 !(group 5-8)

double precision :: real_porosity

! number of wall points
integer :: nWall 

integer, DIMENSION(:), ALLOCATABLE :: vecWall
integer, DIMENSION(:,:), ALLOCATABLE :: dir1, dir2, &
    dir3, dir4, dir5, dir6, dir7, dir8
double precision, DIMENSION(:,:), ALLOCATABLE :: coef1, coef2, coef3, coef4, &
    coef5, coef6, coef7, coef8
double precision, DIMENSION(:,:), ALLOCATABLE :: f1w,f2w,f3w,f4w, &
    f5w,f6w,f7w,f8w

contains 
    ! should be called after calling MPIParams::setupVirtualProcessGrid
    ! such that xl, xlg, xu, xug and yl, ylg, yu, yug has been set already
    subroutine setupPhysicalGrid
        use velocityGrid, only: Nc8, DiffFlux
        implicit none

        ! local vars
        integer :: localid, i, j, k, icount, countVoidP, NneighborFluid
        integer :: bxl, bxu, byl, byu, bzl, bzu ! bound when counting fluid point for sweeping
        integer :: ii, jj, kk
        integer :: nCorner
        character(kind=1) :: ctemp

        ds = 1.d0/(Nx-1-2*fluidlayer) 

        ! set the extend and sizes
        Nxtotal = xug - xlg + 1
        Nytotal = yug - ylg + 1
        Nztotal = zug - zlg + 1
        Nxytotal = Nxtotal*Nytotal
        Nxsub = xu - xl + 1
        Nysub = yu - yl + 1
        Nzsub = zu - zl + 1
        Ntotal = Nxtotal * Nytotal * Nztotal

        allocate(array3D(Nx, Ny, Nz)) ! this is global raw geometry data
        !including ghost layer (flag = ghost)
        allocate(array3Dg(xmin-ghostLayers:xmax+ghostLayers, &
                          ymin-ghostLayers:ymax+ghostLayers, &
                          zmin-ghostLayers:zmax+ghostLayers))

        allocate(image(xlg:xug, ylg:yug, zlg:zug)) ! this is local image 
        allocate(which_corner(xlg:xug, ylg:yug, zlg:zug)) 


        ! set global image array
        array3D  = fluid ! set all to fluid first (fluidLayers are fluid)
        array3Dg = fluid ! set to ghost, then set inner ridge

        !-----------------------------------------------------------------
        !read base digital image, NOTE that all raw image dimension are NY^3. 
        !-----------------------------------------------------------------
        Open(IMAGEFILE,file=imageFileName,status='OLD', &
            form='unformatted',ACCESS="STREAM")
        Do k=1,Nz_base
            Do j=1,Ny_base
                Do i=fluidlayer+1,Nx_base-fluidlayer ! may need to be modified 
                    read(IMAGEFILE) ctemp       
                    !l=(i-translate+ghostLayer)+(j-translate+ghostLayer-1)*Nxtotal+(k-translate+ghostLayer-1)*Nxytotal
                    if ((k>=1+translate).AND.(k<=Nz_base+translate) &
                        .AND.(j>=1+translate).AND.(j<=Ny_base+translate) &
                        .AND.(i>=fluidlayer+1+translate) &
                        .AND.(i<=Nx_base-fluidlayer+translate) ) then
                        if (ichar(ctemp)==0) then
                            array3D(i,j,k) = fluid
                        else if (ichar(ctemp)==1) then
                            array3D(i,j,k) = solid
                        else 
                            write(*,*) "invalid image data file"
                            stop 0
                        endif
                    endif
                Enddo  
            Enddo
        Enddo
        Close(IMAGEFILE)

        !-----------------------------------------------------------------
        !repeat the base block for weak scaling efficiency study. 
        !-----------------------------------------------------------------
        do k = 0, block_repz-1
        do j = 0, block_repy-1
        do i = 0, block_repx-1
            do kk =1, Nz_base
                do jj = 1, Ny_base
                    do ii = 1, Nx_base
                        array3D(i*Nx_base+ii, j*Ny_base+jj, k*Nz_base+kk) &
                            = array3D(ii, jj, kk)
                    enddo
                enddo
            enddo
        enddo
        enddo
        enddo
        
        !-----------------------------------------------------------------   
        ! set ghost and inner of array3Dg
        !-----------------------------------------------------------------
        do k=zmin-ghostLayers,zmax+ghostLayers
            do j=ymin-ghostLayers,ymax+ghostLayers
                do i=xmin-ghostLayers,xmax+ghostLayers
                    if(k<zmin .or. k>zmax .or. j<ymin .or. &
                       j>ymax .or. i<xmin .or. i>xmax) then ! ghost layers
                        array3Dg(i,j,k) = ghost
                    else ! inner region
                        array3Dg(i,j,k) = array3D(i,j,k)
                    endif
                enddo
            enddo
        enddo

        !-----------------------------------------------------------------
        ! drill the small pores (change array3Dg)
        !-----------------------------------------------------------------
        do k=1,Nz
            do j=1,Ny
                do i=2,Nx-1
                    if(array3Dg(i,j,k) == fluid) then
                        if ((array3Dg(i+1,j,k) /=fluid .and. &
                             array3Dg(i-1,j,k) /= fluid) .or. &
                            (array3Dg(i,j+1,k) /=fluid .and. &
                             array3Dg(i,j-1,k) /= fluid) .or. &
                            (array3Dg(i,j,k+1) /=fluid .and. &
                             array3Dg(i,j,k-1) /= fluid) ) then
                            if (j==Ny) then
                                if (k==Nz) then
                                    array3Dg(i:i+1,j-1:j,k-1:k) = fluid
                                else
                                    array3Dg(i:i+1,j-1:j,k:k+1) = fluid
                                endif
                            else if (k==Nz) then
                                array3Dg(i:i+1,j:j+1,k-1:k) = fluid
                            else
                                array3Dg(i:i+1,j:j+1,k:k+1) = fluid
                            endif
                        endif
                    endif !fluid
                enddo !i
            enddo !j
        enddo !k

        !-----------------------------------------------------------------
        ! Remove 1-layer-thickness wall 
        !-----------------------------------------------------------------
        Do k=1,Nz
          Do j=1,Ny
            Do i=1,Nx
              If (array3Dg(i,j,k)==solid) then
                If(   ((array3Dg(i+1,j,k)==fluid).AND.(array3Dg(i-1,j,k)==fluid)) &
                  .OR.((array3Dg(i,j+1,k)==fluid).AND.(array3Dg(i,j-1,k)==fluid))&
                  .OR.((array3Dg(i,j,k+1)==fluid).AND.(array3Dg(i,j,k-1)==fluid))) then
                    array3Dg(i,j,k)=fluid
                Endif       
              Endif
            Enddo
          Enddo
        Enddo  

        !-----------------------------------------------------------------
        ! revoving wall node that near the communication boundary
        ! NOTE: 2017-5-4
        !   the which_corner bug has been found, so no need now
        !-----------------------------------------------------------------
        !                      |                                        
        !       # # # # # # O *|*                                       
        !       # # # # # # O *|*                                       
        !       * # # # # # O *|*                                       
        !       * * O O O O * *|*                                       
        !       * * * * * * * *|*                                       
        !                                                                
        ! NOTE: O type wall points need to be change to fluid
        !       The same role applys to Y/Z commu. boundary
        !-----------------------------------------------------------------
        ! do k=zl,zu
        !   do j=yl,yu
        !     do i=xl,xu
        !       if(array3Dg(i,j,k) == solid) then
        !         if ( (k==(zl+1) .and. array3Dg(i,j,k-1)==fluid) &
        !         .or. (k==(zu-1) .and. array3Dg(i,j,k+1)==fluid) &
        !         .or. (j==(yl+1) .and. array3Dg(i,j-1,k)==fluid) &
        !         .or. (j==(yu-1) .and. array3Dg(i,j+1,k)==fluid) &
        !         .or. (i==(xl+1) .and. array3Dg(i-1,j,k)==fluid) &
        !         .or. (i==(xu-1) .and. array3Dg(i+1,j,k)==fluid) ) then
        !           array3Dg(i,j,k) = fluid
        !         endif
        !       endif
        !     enddo
        !   enddo
        ! enddo


        !-----------------------------------------------------------------
        ! reset array3D, since array3Dg has now been cleaned
        !-----------------------------------------------------------------
        do k = zmin, zmax
            do j = ymin, ymax
                do i = xmin, xmax
                    array3D(i,j,k) = array3Dg(i,j,k)
                enddo
            enddo
        enddo

        ! set local image
        image = fluid
        bxl = xlg
        byl = ylg
        bzl = zlg
        bxu = xug
        byu = yug
        bzu = zug

        if(xl == xmin) bxl = xl !if most west block(processor)
        if(yl == ymin) byl = yl !if most south block
        if(zl == zmin) bzl = zl !if most back  block
        if(xu == xmax) bxu = xu !if most east block(processor)
        if(yu == ymax) byu = yu !if most north block
        if(zu == zmax) bzu = zu !if most front block
        image(bxl:bxu, byl:byu, bzl:bzu) = array3Dg(bxl:bxu, byl:byu, bzl:bzu)

        !Assign the outerboarder layers (called ghost point in serial program)
        If (ghostLayers>0) then
            Do k=zlg,zug
                Do j=ylg,yug
                    Do i=xlg,xug
                        ! test if boundary processor, here the ghost flag doesn't means 
                        ! the communication boundary, but only means the true global 
                        ! outter boarder.
                        If ((i<xmin).OR.(i>xmax).OR.(j<ymin).OR.(j>ymax) &
                            .OR.(k<zmin).OR.(k>zmax)) then 
                            image(i,j,k) = ghost
                        End if
                    Enddo
                Enddo
            Enddo
        End if

        ! set wall points type based on sournding point type(f/s)
        nWall=0 ! count the wall points
        Do k=bzl,bzu
            Do j=byl,byu
                Do i=bxl,bxu
                    !ii = i-xl+1
                    !jj = j-yl+1
                    !localid = (j-ylg)*Nxtotal + i-xlg+1
                    If (array3D(i,j,k)==solid) then
                        NneighborFluid=1 !(1-2-3-4 in D2Q9 corresponding to 2-3-5-7)
                        If (array3Dg(i+1, j, k)==fluid)  NneighborFluid=NneighborFluid*2   !Check neighbor on the East(2)
                        If (array3Dg(i-1, j, k)==fluid)  NneighborFluid=NneighborFluid*3   !Check neighbor on the North(3)
                        If (array3Dg(i, j+1, k)==fluid)  NneighborFluid=NneighborFluid*5   !Check neighbor on the West(5)
                        If (array3Dg(i, j-1, k)==fluid)  NneighborFluid=NneighborFluid*7   !Check neighbor on the South(7)
                        If (array3Dg(i, j, k+1)==fluid)  NneighborFluid=NneighborFluid*9   !Check neighbor on the West(5)
                        If (array3Dg(i, j, k-1)==fluid)  NneighborFluid=NneighborFluid*11   !Check neighbor on the South(7)                        
    
                        SELECT Case (NneighborFluid)
                            CASE (2)
                                image(i,j,k)=WallE
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (3)
                                image(i,j,k)=WallW
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (5)
                                image(i,j,k)=WallN
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (7)
                                image(i,j,k)=WallS
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (9)
                                image(i,j,k)=WallF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (11)
                                image(i,j,k)=WallB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (10)
                                image(i,j,k)=WallEN
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (15)
                                image(i,j,k)=WallWN
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (14)
                                image(i,j,k)=WallES
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (21)
                                image(i,j,k)=WallWS
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (18)
                                image(i,j,k)=WallEF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (22)
                                image(i,j,k)=WallEB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (27)
                                image(i,j,k)=WallWF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (33)
                                image(i,j,k)=WallWB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (45)
                                image(i,j,k)=WallNF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (55)
                                image(i,j,k)=WallNB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (63)
                                image(i,j,k)=WallSF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (77)
                                image(i,j,k)=WallSB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (90)
                                image(i,j,k)=WallENF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (110)
                                image(i,j,k)=WallENB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (135)
                                image(i,j,k)=WallWNF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (165)
                                image(i,j,k)=WallWNB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (126)
                                image(i,j,k)=WallESF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (154)
                                image(i,j,k)=WallESB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (189)
                                image(i,j,k)=WallWSF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            CASE (231)
                                image(i,j,k)=WallWSB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1                                                                                               
                        END SELECT
                    Endif
                Enddo
            Enddo
        Enddo

        PRINT*, "nWall = ", nWall

        ! Create wall-type vectors and record the location of corner points
        ! while the which_corner(i,j,k) array
        ! NOTE: 2017-05-04
        !   Bug found here
        !   Only non-corner wall nodes in ghost layers are not counted
        !
        !vecWall(walli) mark the global id of the walli'th wall in the image
        ALLOCATE(vecWall(nWall))

        nWall=0
        nCorner=0
        do k=zlg, zug
            do j=ylg, yug
                do i=xlg, xug
                    localid = (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                    if(any(allWallLabels(:) == image(i,j,k))) then
                        if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                           k .ge. zl .and. k .le. zu) then
                            nWall=nWall+1
                            vecWall(nWall) = localid
                        endif
                        if(any(onlyCornerLabels(:) == image(i,j,k))) then
                            nCorner=nCorner+1
                            which_corner(i,j,k)=nCorner
                        endif
                    endif
                enddo
            enddo
        enddo                           

        ALLOCATE(f1w(nCorner, Nc8),f2w(nCorner, Nc8),&
                 f3w(nCorner, Nc8),f4w(nCorner, Nc8))
        ALLOCATE(f5w(nCorner, Nc8),f6w(nCorner, Nc8),&
                 f7w(nCorner, Nc8),f8w(nCorner, Nc8))


        !Direction 1
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xl == xmin) bxl = xl+1 !if most west block(processor)
        if(yl == ymin) byl = yl+1 !if most south block
        if(zl == zmin) bzl = zl+1!if most back  block
        Nstencil1=0
        Do k=bzl,bzu
            Do j=byl,byu
                Do i=bxl,bxu
                    If (image(i,j,k)==fluid) then
                        Nstencil1=Nstencil1+1
                    End if
                End do
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir1(3, Nstencil1))
        icount=0
        Do k=bzl,bzu
            Do j=byl,byu
                Do i=bxl,bxu
                    If (image(i,j,k)==fluid) then
                        icount=icount+1
                        dir1(:,icount)=(/i,j,k/)
                    End if
                End do
            End do
        End do

        !Direction 2
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xu == xmax) bxu = xu-1 !if most west block(processor)
        if(yl == ymin) byl = yl+1 !if most south block
        if(zl == zmin) bzl = zl+1!if most back  block
        Nstencil2=0
        Do k=bzl,bzu
            Do j=byl,byu
                Do i=bxu,bxl,-1
                    If (image(i,j,k)==fluid) then
                        Nstencil2=Nstencil2+1
                    End if
                End do
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir2(3, Nstencil2))
        icount=0
        Do k=bzl,bzu
            Do j=byl,byu
                Do i=bxu,bxl,-1
                    If (image(i,j,k)==fluid) then
                        icount=icount+1
                        dir2(:,icount)=(/i,j,k/)
                    End if
                End do
            End do
        End do

        !Direction 3
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xu == xmax) bxu = xu-1 !if most west block(processor)
        if(yu == ymax) byu = yu-1 !if most south block
        if(zl == zmin) bzl = zl+1!if most back  block
        Nstencil3=0
        Do k=bzl,bzu
            Do j=byu,byl,-1
                Do i=bxu,bxl,-1
                    If (image(i,j,k)==fluid) then
                        Nstencil3=Nstencil3+1
                    End if
                End do
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir3(3, Nstencil3))
        icount=0
        Do k=bzl,bzu
            Do j=byu,byl,-1
                Do i=bxu,bxl,-1
                    If (image(i,j,k)==fluid) then
                        icount=icount+1
                        dir3(:,icount)=(/i,j,k/)
                    End if
                End do
            End do
        End do

        !Direction 4
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xl == xmin) bxl = xl+1 !if most west block(processor)
        if(yu == ymax) byu = yu-1 !if most south block
        if(zl == zmin) bzl = zl+1 !if most back  block
        Nstencil4=0
        Do k=bzl,bzu
            Do j=byu,byl,-1
                Do i=bxl,bxu
                    If (image(i,j,k)==fluid) then
                        Nstencil4=Nstencil4+1
                    End if
                End do
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir4(3, Nstencil4))
        icount=0
        Do k=bzl,bzu
            Do j=byu,byl,-1
                Do i=bxl,bxu
                    If (image(i,j,k)==fluid) then
                        icount=icount+1
                        dir4(:,icount)=(/i,j,k/)
                    End if
                End do
            End do
        End do

        !Direction 5
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xl == xmin) bxl = xl+1 !if most west block(processor)
        if(yl == ymin) byl = yl+1 !if most south block
        if(zu == zmax) bzu = zu-1!if most back  block
        Nstencil5=0
        Do k=bzu,bzl,-1
            Do j=byl,byu
                Do i=bxl,bxu
                    If (image(i,j,k)==fluid) then
                        Nstencil5=Nstencil5+1
                    End if
                End do
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir5(3, Nstencil5))
        icount=0
        Do k=bzu,bzl,-1
            Do j=byl,byu
                Do i=bxl,bxu
                    If (image(i,j,k)==fluid) then
                        icount=icount+1
                        dir5(:,icount)=(/i,j,k/)
                    End if
                End do
            End do
        End do

        !Direction 6
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xu == xmax) bxu = xu-1 !if most west block(processor)
        if(yl == ymin) byl = yl+1 !if most south block
        if(zu == zmax) bzu = zu-1!if most back  block
        Nstencil6=0
        Do k=bzu,bzl,-1
            Do j=byl,byu
                Do i=bxu,bxl,-1
                    If (image(i,j,k)==fluid) then
                        Nstencil6=Nstencil6+1
                    End if
                End do
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir6(3, Nstencil6))
        icount=0
        Do k=bzu,bzl,-1
            Do j=byl,byu
                Do i=bxu,bxl,-1
                    If (image(i,j,k)==fluid) then
                        icount=icount+1
                        dir6(:,icount)=(/i,j,k/)
                    End if
                End do
            End do
        End do

        !Direction 7
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xu == xmax) bxu = xu-1 !if most west block(processor)
        if(yu == ymax) byu = yu-1 !if most south block
        if(zu == zmax) bzu = zu-1!if most back  block
        Nstencil7=0
        Do k=bzu,bzl,-1
            Do j=byu,byl,-1
                Do i=bxu,bxl,-1
                    If (image(i,j,k)==fluid) then
                        Nstencil7=Nstencil7+1
                    End if
                End do
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir7(3, Nstencil7))
        icount=0
        Do k=bzu,bzl,-1
            Do j=byu,byl,-1
                Do i=bxu,bxl,-1
                    If (image(i,j,k)==fluid) then
                        icount=icount+1
                        dir7(:,icount)=(/i,j,k/)
                    End if
                End do
            End do
        End do

        !Direction 8
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xl == xmin) bxl = xl+1
        if(yu == ymax) byu = yu-1
        if(zu == zmax) bzu = zu-1
        Nstencil8=0
        Do k=bzu,bzl,-1
            Do j=byu,byl,-1
                Do i=bxl,bxu
                    If (image(i,j,k)==fluid) then
                        Nstencil8=Nstencil8+1
                    End if
                End do
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir8(3, Nstencil8))
        icount=0
        Do k=bzu,bzl,-1
            Do j=byu,byl,-1
                Do i=bxl,bxu
                    If (image(i,j,k)==fluid) then
                        icount=icount+1
                        dir8(:,icount)=(/i,j,k/)
                    End if
                End do
            End do
        End do

        !construct the deferencial coefficients
        ALLOCATE(coef1(Nstencil1,9),coef2(Nstencil2,9),&
            coef3(Nstencil3,9),coef4(Nstencil4,9))
        ALLOCATE(coef5(Nstencil5,9),coef6(Nstencil6,9),&
            coef7(Nstencil7,9),coef8(Nstencil8,9))
        
        Do i=1,Nstencil1
            ii=dir1(1,i) 
            jj=dir1(2,i)
            kk=dir1(3,i)
            !2nd order of accuracy
            coef1(i,1)= 1.5d0/ds !x0n
            coef1(i,2)= 2.d0/ds  !x1n
            coef1(i,3)=-0.5d0/ds   !x2n
            coef1(i,4)= 1.5d0/ds !y0n
            coef1(i,5)= 2.d0/ds  !y1n
            coef1(i,6)=-0.5d0/ds   !y2n
            coef1(i,7)= 1.5d0/ds !z0n
            coef1(i,8)= 2.d0/ds  !z1n
            coef1(i,9)=-0.5d0/ds   !z2n
            
            !1st order of accuracy in x
            if ((image(ii-2,jj,kk)==ghost)  .OR.(image(ii-1,jj,kk)==wallE)  &
            .OR.(image(ii-1,jj,kk)==wallES) .OR.(image(ii-1,jj,kk)==wallEN) &
            .OR.(image(ii-1,jj,kk)==wallEF) .OR.(image(ii-1,jj,kk)==wallEB) &
            .OR.(image(ii-1,jj,kk)==wallENF).OR.(image(ii-1,jj,kk)==wallENB)&
            .OR.(image(ii-1,jj,kk)==wallESF).OR.(image(ii-1,jj,kk)==wallESB)) then
                coef1(i,1)=1.d0/ds  !x0n
                coef1(i,2)=1.d0/ds  !x1n
                coef1(i,3)=0.d0     !x2n
            end if
            !1st order of accuracy in y
            if ((image(ii,jj-2,kk)==ghost)  .OR.(image(ii,jj-1,kk)==wallN)  &
            .OR.(image(ii,jj-1,kk)==wallEN) .OR.(image(ii,jj-1,kk)==wallWN) &
            .OR.(image(ii,jj-1,kk)==wallNF) .OR.(image(ii,jj-1,kk)==wallNB) &
            .OR.(image(ii,jj-1,kk)==wallENF).OR.(image(ii,jj-1,kk)==wallENB)&
            .OR.(image(ii,jj-1,kk)==wallWNF).OR.(image(ii,jj-1,kk)==wallWNB)) then
                coef1(i,4)=1.d0/ds  !y0n
                coef1(i,5)=1.d0/ds  !y1n
                coef1(i,6)=0.d0     !y2n
            end if
            !1st order of accuracy in z
            if ((image(ii,jj,kk-2)==ghost)  .OR.(image(ii,jj,kk-1)==wallF)  &
            .OR.(image(ii,jj,kk-1)==wallNF) .OR.(image(ii,jj,kk-1)==wallSF) &
            .OR.(image(ii,jj,kk-1)==wallEF) .OR.(image(ii,jj,kk-1)==wallWF) &
            .OR.(image(ii,jj,kk-1)==wallENF).OR.(image(ii,jj,kk-1)==wallWNF)&
            .OR.(image(ii,jj,kk-1)==wallESF).OR.(image(ii,jj,kk-1)==wallWSF)) then
                coef1(i,7)=1.d0/ds  !z0n
                coef1(i,8)=1.d0/ds  !z1n
                coef1(i,9)=0.d0     !z2n
            end if
        End do ! end of ceof for dir 1

        Do i=1,Nstencil2
            ii=dir2(1,i) 
            jj=dir2(2,i)
            kk=dir2(3,i)
            !2nd order of accuracy
            coef2(i,1)=-1.5d0/ds !x0n
            coef2(i,2)=-2.d0/ds  !x1n
            coef2(i,3)=0.5d0/ds   !x2n
            coef2(i,4)=1.5d0/ds !y0n
            coef2(i,5)=2.d0/ds  !y1n
            coef2(i,6)=-0.5d0/ds   !y2n
            coef2(i,7)=1.5d0/ds !z0n
            coef2(i,8)=2.d0/ds  !z1n
            coef2(i,9)=-0.5d0/ds   !z2n
            !1st order of accuracy in x
            if ((image(ii+2,jj,kk)==ghost)  .OR.(image(ii+1,jj,kk)==wallW)  &
            .OR.(image(ii+1,jj,kk)==wallWS) .OR.(image(ii+1,jj,kk)==wallWN) &
            .OR.(image(ii+1,jj,kk)==wallWF) .OR.(image(ii+1,jj,kk)==wallWB) &
            .OR.(image(ii+1,jj,kk)==wallWNF).OR.(image(ii+1,jj,kk)==wallWNB)&
            .OR.(image(ii+1,jj,kk)==wallWSF).OR.(image(ii+1,jj,kk)==wallWSB)) then
                coef2(i,1)=-1.d0/ds  !x0n
                coef2(i,2)=-1.d0/ds  !x1n
                coef2(i,3)=0.d0     !x2n
            end if
            !1st order of accuracy in y
            if ((image(ii,jj-2,kk)==ghost)  .OR.(image(ii,jj-1,kk)==wallN)  &
            .OR.(image(ii,jj-1,kk)==wallEN) .OR.(image(ii,jj-1,kk)==wallWN) &
            .OR.(image(ii,jj-1,kk)==wallNF) .OR.(image(ii,jj-1,kk)==wallNB) &
            .OR.(image(ii,jj-1,kk)==wallENF).OR.(image(ii,jj-1,kk)==wallENB)&
            .OR.(image(ii,jj-1,kk)==wallWNF).OR.(image(ii,jj-1,kk)==wallWNB)) then
                coef2(i,4)=1.d0/ds  !y0n
                coef2(i,5)=1.d0/ds  !y1n
                coef2(i,6)=0.d0     !y2n
            end if
            !1st order of accuracy in z
            if ((image(ii,jj,kk-2)==ghost)  .OR.(image(ii,jj,kk-1)==wallF)  &
            .OR.(image(ii,jj,kk-1)==wallNF) .OR.(image(ii,jj,kk-1)==wallSF) &
            .OR.(image(ii,jj,kk-1)==wallEF) .OR.(image(ii,jj,kk-1)==wallWF) &
            .OR.(image(ii,jj,kk-1)==wallENF).OR.(image(ii,jj,kk-1)==wallWNF)&
            .OR.(image(ii,jj,kk-1)==wallESF).OR.(image(ii,jj,kk-1)==wallWSF)) then
                coef2(i,7)=1.d0/ds  !z0n
                coef2(i,8)=1.d0/ds  !z1n
                coef2(i,9)=0.d0     !z2n
            end if
        End do

        Do i=1,Nstencil3
            ii=dir3(1,i) 
            jj=dir3(2,i)
            kk=dir3(3,i)
            !2nd order of accuracy
            coef3(i,1)=-1.5d0/ds !x0n
            coef3(i,2)=-2.d0/ds  !x1n
            coef3(i,3)=0.5d0/ds   !x2n
            coef3(i,4)=-1.5d0/ds !y0n
            coef3(i,5)=-2.d0/ds  !y1n
            coef3(i,6)=0.5d0/ds   !y2n
            coef3(i,7)=1.5d0/ds !z0n
            coef3(i,8)=2.d0/ds  !z1n
            coef3(i,9)=-0.5d0/ds   !z2n
            !1st order of accuracy in x
            if ((image(ii+2,jj,kk)==ghost)  .OR.(image(ii+1,jj,kk)==wallW)  &
            .OR.(image(ii+1,jj,kk)==wallWS) .OR.(image(ii+1,jj,kk)==wallWN) &
            .OR.(image(ii+1,jj,kk)==wallWF) .OR.(image(ii+1,jj,kk)==wallWB) &
            .OR.(image(ii+1,jj,kk)==wallWNF).OR.(image(ii+1,jj,kk)==wallWNB)&
            .OR.(image(ii+1,jj,kk)==wallWSF).OR.(image(ii+1,jj,kk)==wallWSB)) then
                coef3(i,1)=-1.d0/ds  !x0n
                coef3(i,2)=-1.d0/ds  !x1n
                coef3(i,3)=0.d0     !x2n
            end if
            !1st order of accuracy in y
            if ((image(ii,jj+2,kk)==ghost)  .OR.(image(ii,jj+1,kk)==wallS)  &
            .OR.(image(ii,jj+1,kk)==wallES) .OR.(image(ii,jj+1,kk)==wallWS) &
            .OR.(image(ii,jj+1,kk)==wallSF) .OR.(image(ii,jj+1,kk)==wallSB) &
            .OR.(image(ii,jj+1,kk)==wallESF).OR.(image(ii,jj+1,kk)==wallESB)&
            .OR.(image(ii,jj+1,kk)==wallWSF).OR.(image(ii,jj+1,kk)==wallWSB)) then
                coef3(i,4)=-1.d0/ds  !y0n
                coef3(i,5)=-1.d0/ds  !y1n
                coef3(i,6)=0.d0     !y2n
            end if
            !1st order of accuracy in z
            if ((image(ii,jj,kk-2)==ghost)  .OR.(image(ii,jj,kk-1)==wallF)  &
            .OR.(image(ii,jj,kk-1)==wallNF) .OR.(image(ii,jj,kk-1)==wallSF) &
            .OR.(image(ii,jj,kk-1)==wallEF) .OR.(image(ii,jj,kk-1)==wallWF) &
            .OR.(image(ii,jj,kk-1)==wallENF).OR.(image(ii,jj,kk-1)==wallWNF)&
            .OR.(image(ii,jj,kk-1)==wallESF).OR.(image(ii,jj,kk-1)==wallWSF)) then
                coef3(i,7)=1.d0/ds  !z0n
                coef3(i,8)=1.d0/ds  !z1n
                coef3(i,9)=0.d0     !z2n
            end if
        End do
        
        Do i=1,Nstencil4
            ii=dir4(1,i) 
            jj=dir4(2,i)
            kk=dir4(3,i)
            !2nd order of accuracy
            coef4(i,1)=1.5d0/ds !x0n
            coef4(i,2)=2.d0/ds  !x1n
            coef4(i,3)=-0.5d0/ds   !x2n
            coef4(i,4)=-1.5d0/ds !y0n
            coef4(i,5)=-2.d0/ds  !y1n
            coef4(i,6)=0.5d0/ds   !y2n
            coef4(i,7)=1.5d0/ds !z0n
            coef4(i,8)=2.d0/ds  !z1n
            coef4(i,9)=-0.5d0/ds   !z2n
            !1st order of accuracy in x
            if ((image(ii-2,jj,kk)==ghost)  .OR.(image(ii-1,jj,kk)==wallE)   &
            .OR.(image(ii-1,jj,kk)==wallES) .OR.(image(ii-1,jj,kk)==wallEN)  &
            .OR.(image(ii-1,jj,kk)==wallEF) .OR.(image(ii-1,jj,kk)==wallEB)  &
            .OR.(image(ii-1,jj,kk)==wallENF).OR.(image(ii-1,jj,kk)==wallENB) &
            .OR.(image(ii-1,jj,kk)==wallESF).OR.(image(ii-1,jj,kk)==wallESB))  then
                coef4(i,1)=1.d0/ds  !x0n
                coef4(i,2)=1.d0/ds  !x1n
                coef4(i,3)=0.d0     !x2n
            end if
            !1st order of accuracy in y
            if ((image(ii,jj+2,kk)==ghost)  .OR.(image(ii,jj+1,kk)==wallS)  &
            .OR.(image(ii,jj+1,kk)==wallES) .OR.(image(ii,jj+1,kk)==wallWS) &
            .OR.(image(ii,jj+1,kk)==wallSF) .OR.(image(ii,jj+1,kk)==wallSB) &
            .OR.(image(ii,jj+1,kk)==wallESF).OR.(image(ii,jj+1,kk)==wallESB)&
            .OR.(image(ii,jj+1,kk)==wallWSF).OR.(image(ii,jj+1,kk)==wallWSB)) then
                coef4(i,4)=-1.d0/ds  !y0n
                coef4(i,5)=-1.d0/ds  !y1n
                coef4(i,6)=0.d0     !y2n
            end if
            !1st order of accuracy in z
            if ((image(ii,jj,kk-2)==ghost)  .OR.(image(ii,jj,kk-1)==wallF)  &
            .OR.(image(ii,jj,kk-1)==wallNF) .OR.(image(ii,jj,kk-1)==wallSF) &
            .OR.(image(ii,jj,kk-1)==wallEF) .OR.(image(ii,jj,kk-1)==wallWF) &
            .OR.(image(ii,jj,kk-1)==wallENF).OR.(image(ii,jj,kk-1)==wallWNF)&
            .OR.(image(ii,jj,kk-1)==wallESF).OR.(image(ii,jj,kk-1)==wallWSF)) then
                coef4(i,7)=1.d0/ds  !z0n
                coef4(i,8)=1.d0/ds  !z1n
                coef4(i,9)=0.d0     !z2n
            end if
        End do
        
        Do i=1,Nstencil5
            ii=dir5(1,i) 
            jj=dir5(2,i)
            kk=dir5(3,i)
            !2nd order of accuracy
            coef5(i,1)=1.5d0/ds !x0n
            coef5(i,2)=2.d0/ds  !x1n
            coef5(i,3)=-0.5d0/ds   !x2n
            coef5(i,4)=1.5d0/ds !y0n
            coef5(i,5)=2.d0/ds  !y1n
            coef5(i,6)=-0.5d0/ds   !y2n
            coef5(i,7)=-1.5d0/ds !z0n
            coef5(i,8)=-2.d0/ds  !z1n
            coef5(i,9)=0.5d0/ds   !z2n
            !1st order of accuracy in x
            if ((image(ii-2,jj,kk)==ghost)  .OR.(image(ii-1,jj,kk)==wallE)  &
            .OR.(image(ii-1,jj,kk)==wallES) .OR.(image(ii-1,jj,kk)==wallEN) &
            .OR.(image(ii-1,jj,kk)==wallEF) .OR.(image(ii-1,jj,kk)==wallEB) &
            .OR.(image(ii-1,jj,kk)==wallENF).OR.(image(ii-1,jj,kk)==wallENB)&
            .OR.(image(ii-1,jj,kk)==wallESF).OR.(image(ii-1,jj,kk)==wallESB)) then
                coef5(i,1)=1.d0/ds  !x0n
                coef5(i,2)=1.d0/ds  !x1n
                coef5(i,3)=0.d0     !x2n
            end if
            !1st order of accuracy in y
            if ((image(ii,jj-2,kk)==ghost)  .OR.(image(ii,jj-1,kk)==wallN)  &
            .OR.(image(ii,jj-1,kk)==wallEN) .OR.(image(ii,jj-1,kk)==wallWN) &
            .OR.(image(ii,jj-1,kk)==wallNF) .OR.(image(ii,jj-1,kk)==wallNB) &
            .OR.(image(ii,jj-1,kk)==wallENF).OR.(image(ii,jj-1,kk)==wallENB)&
            .OR.(image(ii,jj-1,kk)==wallWNF).OR.(image(ii,jj-1,kk)==wallWNB)) then
                coef5(i,4)=1.d0/ds  !y0n
                coef5(i,5)=1.d0/ds  !y1n
                coef5(i,6)=0.d0     !y2n
            end if
            !1st order of accuracy in z
            if ((image(ii,jj,kk+2)==ghost)  .OR.(image(ii,jj,kk+1)==wallB)  &
            .OR.(image(ii,jj,kk+1)==wallNB) .OR.(image(ii,jj,kk+1)==wallSB) &
            .OR.(image(ii,jj,kk+1)==wallEB) .OR.(image(ii,jj,kk+1)==wallWB) &
            .OR.(image(ii,jj,kk+1)==wallENB).OR.(image(ii,jj,kk+1)==wallWNB)&
            .OR.(image(ii,jj,kk+1)==wallESB).OR.(image(ii,jj,kk+1)==wallWSB)) then
                coef5(i,7)=-1.d0/ds  !z0n
                coef5(i,8)=-1.d0/ds  !z1n
                coef5(i,9)=0.d0     !z2n
            end if
        End do

        Do i=1,Nstencil6
            ii=dir6(1,i) 
            jj=dir6(2,i)
            kk=dir6(3,i)
            !2nd order of accuracy
            coef6(i,1)=-1.5d0/ds !x0n
            coef6(i,2)=-2.d0/ds  !x1n
            coef6(i,3)=0.5d0/ds   !x2n
            coef6(i,4)=1.5d0/ds !y0n
            coef6(i,5)=2.d0/ds  !y1n
            coef6(i,6)=-0.5d0/ds   !y2n
            coef6(i,7)=-1.5d0/ds !z0n
            coef6(i,8)=-2.d0/ds  !z1n
            coef6(i,9)=0.5d0/ds   !z2n
            !1st order of accuracy in x
            if ((image(ii+2,jj,kk)==ghost)  .OR.(image(ii+1,jj,kk)==wallW)  &
            .OR.(image(ii+1,jj,kk)==wallWS) .OR.(image(ii+1,jj,kk)==wallWN) &
            .OR.(image(ii+1,jj,kk)==wallWF) .OR.(image(ii+1,jj,kk)==wallWB) &
            .OR.(image(ii+1,jj,kk)==wallWNF).OR.(image(ii+1,jj,kk)==wallWNB)&
            .OR.(image(ii+1,jj,kk)==wallWSF).OR.(image(ii+1,jj,kk)==wallWSB)) then
                coef6(i,1)=-1.d0/ds  !x0n
                coef6(i,2)=-1.d0/ds  !x1n
                coef6(i,3)=0.d0     !x2n
            end if
            !1st order of accuracy in y
            if ((image(ii,jj-2,kk)==ghost)  .OR.(image(ii,jj-1,kk)==wallN)  &
            .OR.(image(ii,jj-1,kk)==wallEN) .OR.(image(ii,jj-1,kk)==wallWN) &
            .OR.(image(ii,jj-1,kk)==wallNF) .OR.(image(ii,jj-1,kk)==wallNB) &
            .OR.(image(ii,jj-1,kk)==wallENF).OR.(image(ii,jj-1,kk)==wallENB)&
            .OR.(image(ii,jj-1,kk)==wallWNF).OR.(image(ii,jj-1,kk)==wallWNB))  then
                coef6(i,4)=1.d0/ds  !y0n
                coef6(i,5)=1.d0/ds  !y1n
                coef6(i,6)=0.d0     !y2n
            end if
            !1st order of accuracy in z
            if ((image(ii,jj,kk+2)==ghost)  .OR.(image(ii,jj,kk+1)==wallB)  &
            .OR.(image(ii,jj,kk+1)==wallNB) .OR.(image(ii,jj,kk+1)==wallSB) &
            .OR.(image(ii,jj,kk+1)==wallEB) .OR.(image(ii,jj,kk+1)==wallWB) &
            .OR.(image(ii,jj,kk+1)==wallENB).OR.(image(ii,jj,kk+1)==wallWNB)&
            .OR.(image(ii,jj,kk+1)==wallESB).OR.(image(ii,jj,kk+1)==wallWSB)) then
                coef6(i,7)=-1.d0/ds  !z0n
                coef6(i,8)=-1.d0/ds  !z1n
                coef6(i,9)=0.d0     !z2n
            end if
        End do

        Do i=1,Nstencil7
            ii=dir7(1,i) 
            jj=dir7(2,i)
            kk=dir7(3,i)
            !2nd order of accuracy
            coef7(i,1)=-1.5d0/ds !x0n
            coef7(i,2)=-2.d0/ds  !x1n
            coef7(i,3)=0.5d0/ds   !x2n
            coef7(i,4)=-1.5d0/ds !y0n
            coef7(i,5)=-2.d0/ds  !y1n
            coef7(i,6)=0.5d0/ds   !y2
            coef7(i,7)=-1.5d0/ds !z0n
            coef7(i,8)=-2.d0/ds  !z1n
            coef7(i,9)=0.5d0/ds   !z2n
            !1st order of accuracy in x
            if ((image(ii+2,jj,kk)==ghost)  .OR.(image(ii+1,jj,kk)==wallW)  &
            .OR.(image(ii+1,jj,kk)==wallWS) .OR.(image(ii+1,jj,kk)==wallWN) &
            .OR.(image(ii+1,jj,kk)==wallWF) .OR.(image(ii+1,jj,kk)==wallWB) &
            .OR.(image(ii+1,jj,kk)==wallWNF).OR.(image(ii+1,jj,kk)==wallWNB)&
            .OR.(image(ii+1,jj,kk)==wallWSF).OR.(image(ii+1,jj,kk)==wallWSB)) then
                coef7(i,1)=-1.d0/ds  !x0n
                coef7(i,2)=-1.d0/ds  !x1n
                coef7(i,3)=0.d0     !x2n
            end if
            !1st order of accuracy in y
            if ((image(ii,jj+2,kk)==ghost)  .OR.(image(ii,jj+1,kk)==wallS)  &
            .OR.(image(ii,jj+1,kk)==wallES) .OR.(image(ii,jj+1,kk)==wallWS) &
            .OR.(image(ii,jj+1,kk)==wallSF) .OR.(image(ii,jj+1,kk)==wallSB) &
            .OR.(image(ii,jj+1,kk)==wallESF).OR.(image(ii,jj+1,kk)==wallESB)&
            .OR.(image(ii,jj+1,kk)==wallWSF).OR.(image(ii,jj+1,kk)==wallWSB)) then
                coef7(i,4)=-1.d0/ds  !y0n
                coef7(i,5)=-1.d0/ds  !y1n
                coef7(i,6)=0.d0     !y2n
            end if
            !1st order of accuracy in z
            if ((image(ii,jj,kk+2)==ghost)  .OR.(image(ii,jj,kk+1)==wallB)  &
            .OR.(image(ii,jj,kk+1)==wallNB) .OR.(image(ii,jj,kk+1)==wallSB) &
            .OR.(image(ii,jj,kk+1)==wallEB) .OR.(image(ii,jj,kk+1)==wallWB) &
            .OR.(image(ii,jj,kk+1)==wallENB).OR.(image(ii,jj,kk+1)==wallWNB)&
            .OR.(image(ii,jj,kk+1)==wallESB).OR.(image(ii,jj,kk+1)==wallWSB)) then
                coef7(i,7)=-1.d0/ds  !z0n
                coef7(i,8)=-1.d0/ds  !z1n
                coef7(i,9)=0.d0     !z2n
            end if
        End do


        Do i=1,Nstencil8
            ii=dir8(1,i) 
            jj=dir8(2,i)
            kk=dir8(3,i)
            !2nd order of accuracy
            coef8(i,1)=1.5d0/ds !x0n
            coef8(i,2)=2.d0/ds  !x1n
            coef8(i,3)=-0.5d0/ds   !x2n
            coef8(i,4)=-1.5d0/ds !y0n
            coef8(i,5)=-2.d0/ds  !y1n
            coef8(i,6)=0.5d0/ds   !y2n
            coef8(i,7)=-1.5d0/ds !z0n
            coef8(i,8)=-2.d0/ds  !z1n
            coef8(i,9)=0.5d0/ds   !z2n
            !1st order of accuracy in x
            if ((image(ii-2,jj,kk)==ghost)  .OR.(image(ii-1,jj,kk)==wallE)  &
            .OR.(image(ii-1,jj,kk)==wallES) .OR.(image(ii-1,jj,kk)==wallEN) &
            .OR.(image(ii-1,jj,kk)==wallEF) .OR.(image(ii-1,jj,kk)==wallEB) &
            .OR.(image(ii-1,jj,kk)==wallENF).OR.(image(ii-1,jj,kk)==wallENB)&
            .OR.(image(ii-1,jj,kk)==wallESF).OR.(image(ii-1,jj,kk)==wallESB)) then
                coef8(i,1)=1.d0/ds  !x0n
                coef8(i,2)=1.d0/ds  !x1n
                coef8(i,3)=0.d0     !x2n
            end if
            !1st order of accuracy in y
            if ((image(ii,jj+2,kk)==ghost)  .OR.(image(ii,jj+1,kk)==wallS)  &
            .OR.(image(ii,jj+1,kk)==wallES) .OR.(image(ii,jj+1,kk)==wallWS) &
            .OR.(image(ii,jj+1,kk)==wallSF) .OR.(image(ii,jj+1,kk)==wallSB) &
            .OR.(image(ii,jj+1,kk)==wallESF).OR.(image(ii,jj+1,kk)==wallESB)&
            .OR.(image(ii,jj+1,kk)==wallWSF).OR.(image(ii,jj+1,kk)==wallWSB)) then
                coef8(i,4)=-1.d0/ds  !y0n
                coef8(i,5)=-1.d0/ds  !y1n
                coef8(i,6)=0.d0     !y2n
            end if
            !1st order of accuracy in z
            if ((image(ii,jj,kk+2)==ghost  ).OR.(image(ii,jj,kk+1)==wallB)  &
            .OR.(image(ii,jj,kk+1)==wallNB ).OR.(image(ii,jj,kk+1)==wallSB) &
            .OR.(image(ii,jj,kk+1)==wallEB ).OR.(image(ii,jj,kk+1)==wallWB) &
            .OR.(image(ii,jj,kk+1)==wallENB).OR.(image(ii,jj,kk+1)==wallWNB)&
            .OR.(image(ii,jj,kk+1)==wallESB).OR.(image(ii,jj,kk+1)==wallWSB)) then
                coef8(i,7)=-1.d0/ds  !z0n
                coef8(i,8)=-1.d0/ds  !z1n
                coef8(i,9)=0.d0     !z2n
            end if
        End do

        f1w=0.d0
        f2w=0.d0
        f3w=0.d0
        f4w=0.d0
        f5w=0.d0
        f6w=0.d0
        f7w=0.d0
        f8w=0.d0
    end subroutine setupPhysicalGrid
end module physicalGrid
