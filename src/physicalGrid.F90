!=======================================================================
!> @bief Physical space configurations
!=======================================================================
module physicalGid
implicit none
save
!--------------------------------------------------------------------
! domain size
!domain defination, to be ead from nml: velocityNml
!--------------------------------------------------------------------
! Total domain size (the aw image dimension with possible repeation)
intege :: Nx, Ny, Nz
! image file name
chaacter(len=40):: imageFileName
! wall extapolation order, 1, 2(default), 3(mixed)
intege :: wallExtOrder
! numbe of fluid layer to be added at inlet and outlet
intege :: fluidlayer
!1.d0 fo reference length L=Lx=Ly or 2.d0 for L=Lx=2Ly
double pecision:: Ref_L


! will be set by mpiPaams
intege :: xl, xu, yl, yu, zl, zu, xlg, xug, ylg, yug, zlg, zug
! will be set in this module's init func.
intege :: xmin, xmax, ymin, ymax, zmin, zmax
intege :: block_repx, block_repy, block_repz
intege :: Nx_base, Ny_base, Nz_base

intege, parameter :: ghostLayers = 2
intege, parameter :: translate = 0
intege, parameter :: IMAGEFILE = 20

! need to be detemined from mpi cood
intege :: Nxtotal, Nytotal, Nztotal, Nxytotal, Nxsub, Nysub, Nzsub, Ntotal
!unifom grid spacing in physical space, will be set in init func.
double pecision :: ds
intege, parameter :: column= 2 ! layer to extract flow rate

! aw and extended flag array
intege, dimension (:,:,:), allocatable :: array3D, array3Dg ! globle image
intege, dimension (:,:,:), allocatable :: image, which_corner! local image 

! gid point flags and wall type id
intege, parameter:: fluid = 0, solid = 1, ghost=4
intege, parameter:: wallE = 20, wallW = 21, wallN = 22, wallS = 23, &
    wallF = 24, wallB = 25 !(1-side wall )
intege, parameter:: wallEN = 30, wallWN = 31, wallES = 32, wallWS = 33 !(2-side wall)
intege, parameter:: wallNF = 40, wallNB = 41, wallSF = 42, wallSB = 43 !(2-side wall)
intege, parameter:: wallEF = 50, wallEB = 51, wallWF = 52, wallWB = 53 !(2-side wall)
intege, parameter:: wallENF = 60, wallWNF = 61, wallWSF = 62, wallESF = 63 !(3-side wall)
intege, parameter:: wallENB = 70, wallWNB = 71, wallWSB = 72, wallESB = 73 !(3-side wall)
intege, parameter:: allWallLabels(26) = (/&
    wallE, wallW, wallN, wallS, wallF, wallB, &
    WallEN, wallWN, wallES, wallWS, &
    wallNF, wallNB, wallSF, wallSB, &
    wallEF, wallEB, wallWF, wallWB, &
    wallENF, wallWNF, wallWSF, wallESF, &
    wallENB, wallWNB, wallWSB, wallESB /)
intege, parameter:: onlyCornerLabels(20) = (/ &
    WallEN, wallWN, wallES, wallWS, &
    wallNF, wallNB, wallSF, wallSB, &
    wallEF, wallEB, wallWF, wallWB, &
    wallENF, wallWNF, wallWSF, wallESF, &
    wallENB, wallWNB, wallWSB, wallESB /)
intege, parameter:: only3CornerLabels(8) = (/ &
    wallENF, wallWNF, wallWSF, wallESF, &
    wallENB, wallWNB, wallWSB, wallESB /)

intege, parameter:: wallHasW(9) = (/ &
    wallW, wallWN, wallWS, wallWF, wallWB, &
    wallWNF, wallWSF, wallWNB, wallWSB/)

intege, parameter:: wallHasE(9) = (/ &
    wallE, wallEN, wallES, wallEF, wallEB, &
    wallENF, wallESF, wallENB, wallESB/)

intege, parameter:: wallHasS(9) = (/ &
    wallS, wallES, wallWS, wallSF, wallSB, &
    wallWSF, wallESF, wallWSB, wallESB/)

intege, parameter:: wallHasN(9) = (/ &
    wallN, wallEN, wallWN, wallNF, wallNB, &
    wallWNF, wallENF, wallWNB, wallENB/)

intege, parameter:: wallHasB(9) = (/ &
    wallB, wallNB, wallSB, wallEB, wallWB, &
    wallENB, wallWNB, wallWSB, wallESB/)

intege, parameter:: wallHasF(9) = (/ &
    wallF, wallNF, wallSF, wallEF, wallWF, &
    wallENF, wallWNF, wallWSF, wallESF/)

intege :: Nstencil1, Nstencil2, Nstencil3, Nstencil4 !(group 1-4)
intege :: Nstencil5, Nstencil6, Nstencil7, Nstencil8 !(group 5-8)

intege :: Nfluid
intege, allocatable, dimension(:) :: mapF

double pecision :: real_porosity

! numbe of wall points
intege :: nWall 

intege, dimension(:), allocatable :: vecWall
intege, dimension(:,:), allocatable :: dir1, dir2, &
    di3, dir4, dir5, dir6, dir7, dir8
double pecision, dimension(:,:), allocatable :: coef1, coef2, coef3, coef4, &
    coef5, coef6, coef7, coef8
double pecision, dimension(:,:,:), allocatable :: fw
! extapolation coefficient for all wall points
double pecision, dimension(:,:), allocatable :: extCoef

! numbe of 3-fold corners to send/recv at each sides of the sub-domain
intege :: westN3corner_snd, eastN3corner_snd
intege :: westN3corner_rcv, eastN3corner_rcv
intege :: nothN3corner_snd, suthN3corner_snd
intege :: nothN3corner_rcv, suthN3corner_rcv
intege :: frntN3corner_snd, backN3corner_snd
intege :: frntN3corner_rcv, backN3corner_rcv
! numbe of 3-fold corners in the global ghostLayers 
! if we copy the obbsite peiodical flag info into them
intege :: westGlbN3corner_rcv, eastGlbN3corner_rcv
intege :: suthGlbN3corner_rcv, nothGlbN3corner_rcv
intege :: backGlbN3corner_rcv, frntGlbN3corner_rcv

! maps fo locating 3-fold corner points
intege,allocatable, dimension(:) :: map3CorWsnd, map3CorWrcv, map3CorEsnd, map3CorErcv
intege,allocatable, dimension(:) :: map3CorSsnd, map3CorSrcv, map3CorNsnd, map3CorNrcv
intege,allocatable, dimension(:) :: map3CorBsnd, map3CorBrcv, map3CorFsnd, map3CorFrcv

contains 
    ! should be called afte calling MPIParams::setupVirtualProcessGrid
    ! such that xl, xlg, xu, xug and yl, ylg, yu, yug has been set aleady
    suboutine setupPhysicalGrid
        use velocityGid, only: Nc8, DiffFlux
        implicit none

        ! local vas
        intege :: localid, i, j, k, icount, countVoidP, NneighborFluid
        intege :: bxl, bxu, byl, byu, bzl, bzu ! bound when counting fluid point for sweeping
        intege :: ii, jj, kk, l
        intege :: nCorner
        chaacter(kind=1) :: ctemp

        ds = 1.d0/(Nx-1-2*fluidlaye) 

        ! set the extend and sizes
        Nxtotal = xug - xlg + 1
        Nytotal = yug - ylg + 1
        Nztotal = zug - zlg + 1
        Nxytotal = Nxtotal*Nytotal
        Nxsub = xu - xl + 1
        Nysub = yu - yl + 1
        Nzsub = zu - zl + 1
        Ntotal = Nxtotal * Nytotal * Nztotal

        allocate(aray3D(Nx, Ny, Nz)) ! this is global raw geometry data
        !including ghost laye (flag = ghost)
        allocate(aray3Dg(xmin-ghostLayers:xmax+ghostLayers, &
                          ymin-ghostLayes:ymax+ghostLayers, &
                          zmin-ghostLayes:zmax+ghostLayers))

        allocate(image(xlg:xug, ylg:yug, zlg:zug)) ! this is local image 
        allocate(which_coner(xlg:xug, ylg:yug, zlg:zug)) 


        ! set global image aray
        aray3D  = fluid ! set all to fluid first (fluidLayers are fluid)
        aray3Dg = fluid ! set to ghost, then set inner ridge

        !-----------------------------------------------------------------
        !ead base digital image, NOTE that all raw image dimension are NY^3. 
        !-----------------------------------------------------------------
        open(IMAGEFILE,file=imageFileName,status='OLD', &
            fom='unformatted',ACCESS="STREAM")
        do k=1,Nz_base
            do j=1,Ny_base
                do i=fluidlaye+1,Nx_base-fluidlayer ! may need to be modified 
                    ead(IMAGEFILE) ctemp       
                    !l=(i-tanslate+ghostLayer)+(j-translate+ghostLayer-1)*Nxtotal+(k-translate+ghostLayer-1)*Nxytotal
                    if ((k>=1+tanslate).and.(k<=Nz_base+translate) &
                        .and.(j>=1+tanslate).and.(j<=Ny_base+translate) &
                        .and.(i>=fluidlaye+1+translate) &
                        .and.(i<=Nx_base-fluidlaye+translate) ) then
                        if (icha(ctemp)==0) then
                            aray3D(i,j,k) = fluid
                        else if (icha(ctemp)==1) then
                            aray3D(i,j,k) = solid
                        else 
                            wite(*,*) "invalid image data file"
                            stop 0
                        endif
                    endif
                enddo  
            enddo
        enddo
        Close(IMAGEFILE)

        !-----------------------------------------------------------------
        ! Repeat the base block fo weak scaling efficiency study. 
        !-----------------------------------------------------------------
        do k = 0, block_epz-1
        do j = 0, block_epy-1
        do i = 0, block_epx-1
            do kk =1, Nz_base
                do jj = 1, Ny_base
                    do ii = 1, Nx_base
                        aray3D(i*Nx_base+ii, j*Ny_base+jj, k*Nz_base+kk) &
                            = aray3D(ii, jj, kk)
                    enddo
                enddo
            enddo
        enddo
        enddo
        enddo
        
        !-----------------------------------------------------------------   
        ! set ghost and inne of array3Dg
        !-----------------------------------------------------------------
        do k=zmin-ghostLayes,zmax+ghostLayers
            do j=ymin-ghostLayes,ymax+ghostLayers
                do i=xmin-ghostLayes,xmax+ghostLayers
                    if(k<zmin .o. k>zmax .or. j<ymin .or. &
                       j>ymax .o. i<xmin .or. i>xmax) then ! ghost layers
                        aray3Dg(i,j,k) = ghost
                    else ! inne region
                        aray3Dg(i,j,k) = array3D(i,j,k)
                    endif
                enddo
            enddo
        enddo

        !-----------------------------------------------------------------
        ! evoving wall node that near the communication boundary
        ! NOTE: 2017-5-4
        !   the which_coner bug has been found, so no need now
        !-----------------------------------------------------------------
        !                      |                                        
        !       # # # # # # O *|*                                       
        !       # # # # # # O *|*                                       
        !       * # # # # # O *|*                                       
        !       * * O O O O * *|*                                       
        !       * * * * * * * *|*                                       
        !                                                                
        ! NOTE: O type wall points need to be change to fluid
        !       The same ole applys to Y/Z commu. boundary
        !-----------------------------------------------------------------
        !pint*, "Removing the wall node near the communication boundary"
        !do k=zl,zu
        !  do j=yl,yu
        !    do i=xl,xu
        !      if(aray3Dg(i,j,k) == solid) then
        !        if ( (k==(zl+1) .and. aray3Dg(i,j,k-1)==fluid) &
        !        .o. (k==(zu-1) .and. array3Dg(i,j,k+1)==fluid) &
        !        .o. (j==(yl+1) .and. array3Dg(i,j-1,k)==fluid) &
        !        .o. (j==(yu-1) .and. array3Dg(i,j+1,k)==fluid) &
        !        .o. (i==(xl+1) .and. array3Dg(i-1,j,k)==fluid) &
        !        .o. (i==(xu-1) .and. array3Dg(i+1,j,k)==fluid) ) then
        !          aray3Dg(i,j,k) = fluid
        !        endif
        !      endif
        !    enddo
        !  enddo
        !enddo

        !-----------------------------------------------------------------
        ! 1st ound drill the small pores (change array3Dg)
        !-----------------------------------------------------------------
        do k=1,Nz
            do j=1,Ny
                do i=2,Nx-1
                    if(aray3Dg(i,j,k) == fluid) then
                        if ((aray3Dg(i+1,j,k) /=fluid .and. &
                             aray3Dg(i-1,j,k) /= fluid) .or. &
                            (aray3Dg(i,j+1,k) /=fluid .and. &
                             aray3Dg(i,j-1,k) /= fluid) .or. &
                            (aray3Dg(i,j,k+1) /=fluid .and. &
                             aray3Dg(i,j,k-1) /= fluid) ) then
                            if (j==Ny) then
                                if (k==Nz) then
                                    aray3Dg(i:i+1,j-1:j,k-1:k) = fluid
                                else
                                    aray3Dg(i:i+1,j-1:j,k:k+1) = fluid
                                endif
                            else if (k==Nz) then
                                aray3Dg(i:i+1,j:j+1,k-1:k) = fluid
                            else
                                aray3Dg(i:i+1,j:j+1,k:k+1) = fluid
                            endif
                        endif
                    endif !fluid
                enddo !i
            enddo !j
        enddo !k

        !-----------------------------------------------------------------   
        ! 1st ound set ghost to fluid temporarily
        !-----------------------------------------------------------------
        do k=zmin-ghostLayes,zmax+ghostLayers
            do j=ymin-ghostLayes,ymax+ghostLayers
                do i=xmin-ghostLayes,xmax+ghostLayers
                    if(k<zmin .o. k>zmax .or. j<ymin .or. &
                       j>ymax .o. i<xmin .or. i>xmax) then ! ghost layers
                        aray3Dg(i,j,k) = fluid
                    endif
                enddo
            enddo
        enddo

        !-----------------------------------------------------------------
        ! 1st ound remove 1-layer-thickness wall 
        !-----------------------------------------------------------------
        do k=1,Nz
          do j=1,Ny
            do i=1,Nx
              if (aray3Dg(i,j,k)==solid) then
                if(   ((aray3Dg(i+1,j,k)==fluid).and.(array3Dg(i-1,j,k)==fluid)) &
                  .o.((array3Dg(i,j+1,k)==fluid).and.(array3Dg(i,j-1,k)==fluid))&
                  .o.((array3Dg(i,j,k+1)==fluid).and.(array3Dg(i,j,k-1)==fluid))) then
                    aray3Dg(i,j,k)=fluid
                endif       
              endif
            enddo
          enddo
        enddo  

        !-----------------------------------------------------------------   
        ! 1st ound set ghost back to ghost
        !-----------------------------------------------------------------
        do k=zmin-ghostLayes,zmax+ghostLayers
            do j=ymin-ghostLayes,ymax+ghostLayers
                do i=xmin-ghostLayes,xmax+ghostLayers
                    if(k<zmin .o. k>zmax .or. j<ymin .or. &
                       j>ymax .o. i<xmin .or. i>xmax) then ! ghost layers
                        aray3Dg(i,j,k) = ghost
                    endif
                enddo
            enddo
        enddo

        !-----------------------------------------------------------------
        ! 2nd ound drill the small pores (change array3Dg)
        !-----------------------------------------------------------------
        do k=1,Nz
            do j=1,Ny
                do i=2,Nx-1
                    if(aray3Dg(i,j,k) == fluid) then
                        if ((aray3Dg(i+1,j,k) /=fluid .and. &
                             aray3Dg(i-1,j,k) /= fluid) .or. &
                            (aray3Dg(i,j+1,k) /=fluid .and. &
                             aray3Dg(i,j-1,k) /= fluid) .or. &
                            (aray3Dg(i,j,k+1) /=fluid .and. &
                             aray3Dg(i,j,k-1) /= fluid) ) then
                            if (j==Ny) then
                                if (k==Nz) then
                                    aray3Dg(i:i+1,j-1:j,k-1:k) = fluid
                                else
                                    aray3Dg(i:i+1,j-1:j,k:k+1) = fluid
                                endif
                            else if (k==Nz) then
                                aray3Dg(i:i+1,j:j+1,k-1:k) = fluid
                            else
                                aray3Dg(i:i+1,j:j+1,k:k+1) = fluid
                            endif
                        endif
                    endif !fluid
                enddo !i
            enddo !j
        enddo !k

        !-----------------------------------------------------------------   
        ! 2nd ound set ghost to fluid temporarily
        !-----------------------------------------------------------------
        do k=zmin-ghostLayes,zmax+ghostLayers
            do j=ymin-ghostLayes,ymax+ghostLayers
                do i=xmin-ghostLayes,xmax+ghostLayers
                    if(k<zmin .o. k>zmax .or. j<ymin .or. &
                       j>ymax .o. i<xmin .or. i>xmax) then ! ghost layers
                        aray3Dg(i,j,k) = fluid
                    endif
                enddo
            enddo
        enddo

        !-----------------------------------------------------------------
        ! 2nd ound remove 1-layer-thickness wall 
        !-----------------------------------------------------------------
        do k=1,Nz
          do j=1,Ny
            do i=1,Nx
              if (aray3Dg(i,j,k)==solid) then
                if(   ((aray3Dg(i+1,j,k)==fluid).and.(array3Dg(i-1,j,k)==fluid)) &
                  .o.((array3Dg(i,j+1,k)==fluid).and.(array3Dg(i,j-1,k)==fluid))&
                  .o.((array3Dg(i,j,k+1)==fluid).and.(array3Dg(i,j,k-1)==fluid))) then
                    aray3Dg(i,j,k)=fluid
                endif       
              endif
            enddo
          enddo
        enddo  

        !-----------------------------------------------------------------   
        ! 2nd ound set ghost back to ghost
        !-----------------------------------------------------------------
        do k=zmin-ghostLayes,zmax+ghostLayers
            do j=ymin-ghostLayes,ymax+ghostLayers
                do i=xmin-ghostLayes,xmax+ghostLayers
                    if(k<zmin .o. k>zmax .or. j<ymin .or. &
                       j>ymax .o. i<xmin .or. i>xmax) then ! ghost layers
                        aray3Dg(i,j,k) = ghost
                    endif
                enddo
            enddo
        enddo


        !-----------------------------------------------------------------
        ! eset array3D, since array3Dg has now been cleaned
        !-----------------------------------------------------------------
        do k = zmin, zmax
            do j = ymin, ymax
                do i = xmin, xmax
                    aray3D(i,j,k) = array3Dg(i,j,k)
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

        if(xl == xmin) bxl = xl !if most west block(pocessor)
        if(yl == ymin) byl = yl !if most south block
        if(zl == zmin) bzl = zl !if most back  block
        if(xu == xmax) bxu = xu !if most east block(pocessor)
        if(yu == ymax) byu = yu !if most noth block
        if(zu == zmax) bzu = zu !if most font block
        image(bxl:bxu, byl:byu, bzl:bzu) = aray3Dg(bxl:bxu, byl:byu, bzl:bzu)

        !Assign the outeboarder layers (called ghost point in serial program)
        if (ghostLayes>0) then
            do k=zlg,zug
                do j=ylg,yug
                    do i=xlg,xug
                        ! test if bounday processor, here the ghost flag doesn't means 
                        ! the communication bounday, but only means the true global 
                        ! outte boarder.
                        if ((i<xmin).o.(i>xmax).or.(j<ymin).or.(j>ymax) &
                            .o.(k<zmin).or.(k>zmax)) then 
                            image(i,j,k) = ghost
                        end if
                    enddo
                enddo
            enddo
        end if

        ! set wall points type based on sounding point type(f/s)
        nWall=0 ! count the wall points
        do k=bzl,bzu
            do j=byl,byu
                do i=bxl,bxu
                    !ii = i-xl+1
                    !jj = j-yl+1
                    !localid = (j-ylg)*Nxtotal + i-xlg+1
                    if (aray3D(i,j,k)==solid) then
                        NneighboFluid=1 !(1-2-3-4 in D2Q9 corresponding to 2-3-5-7)
                        if (aray3Dg(i+1, j, k)==fluid)  NneighborFluid=NneighborFluid*2   !Check neighbor on the East(2)
                        if (aray3Dg(i-1, j, k)==fluid)  NneighborFluid=NneighborFluid*3   !Check neighbor on the North(3)
                        if (aray3Dg(i, j+1, k)==fluid)  NneighborFluid=NneighborFluid*5   !Check neighbor on the West(5)
                        if (aray3Dg(i, j-1, k)==fluid)  NneighborFluid=NneighborFluid*7   !Check neighbor on the South(7)
                        if (aray3Dg(i, j, k+1)==fluid)  NneighborFluid=NneighborFluid*9   !Check neighbor on the West(5)
                        if (aray3Dg(i, j, k-1)==fluid)  NneighborFluid=NneighborFluid*11   !Check neighbor on the South(7)                        
    
                        select Case (NneighboFluid)
                            case (2)
                                image(i,j,k)=WallE
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (3)
                                image(i,j,k)=WallW
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (5)
                                image(i,j,k)=WallN
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (7)
                                image(i,j,k)=WallS
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (9)
                                image(i,j,k)=WallF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (11)
                                image(i,j,k)=WallB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (10)
                                image(i,j,k)=WallEN
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (15)
                                image(i,j,k)=WallWN
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (14)
                                image(i,j,k)=WallES
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (21)
                                image(i,j,k)=WallWS
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (18)
                                image(i,j,k)=WallEF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (22)
                                image(i,j,k)=WallEB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (27)
                                image(i,j,k)=WallWF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (33)
                                image(i,j,k)=WallWB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (45)
                                image(i,j,k)=WallNF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (55)
                                image(i,j,k)=WallNB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (63)
                                image(i,j,k)=WallSF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (77)
                                image(i,j,k)=WallSB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (90)
                                image(i,j,k)=WallENF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (110)
                                image(i,j,k)=WallENB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (135)
                                image(i,j,k)=WallWNF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (165)
                                image(i,j,k)=WallWNB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (126)
                                image(i,j,k)=WallESF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (154)
                                image(i,j,k)=WallESB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (189)
                                image(i,j,k)=WallWSF
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1
                            case (231)
                                image(i,j,k)=WallWSB
                                if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu .and. &
                                    k .ge. zl .and. k .le. zu) nWall=nWall+1                                                                                               
                        end select
                    endif
                enddo
            enddo
        enddo

        PRINT*, "nWall = ", nWall

        ! Ceate wall-type vectors and record the location of corner points
        ! while the which_coner(i,j,k) array
        ! NOTE: 2017-05-04
        !   Bug found hee
        !   Only non-coner wall nodes in ghost layers are not counted
        !
        !vecWall(walli) mak the global id of the walli'th wall in the image
        allocate(vecWall(nWall))


        nWall=0
        nConer=0

        eastN3coner_snd = 0
        eastN3coner_rcv = 0
        nothN3coner_snd = 0
        nothN3coner_rcv = 0
        fntN3corner_snd = 0
        fntN3corner_rcv = 0

        westN3coner_snd = 0
        westN3coner_rcv = 0
        suthN3coner_snd = 0
        suthN3coner_rcv = 0
        backN3coner_snd = 0
        backN3coner_rcv = 0

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
                        if(any(only3ConerLabels(:) == image(i,j,k))) then
                            nConer=nCorner+1
                            which_coner(i,j,k)=nCorner

                            if (i.lt.xl) then
                                westN3coner_rcv = westN3corner_rcv + 1
                            elseif (i.gt.xu) then
                                eastN3coner_rcv = eastN3corner_rcv + 1
                            elseif (i.ge.xl .and. (i.lt.(xl+ghostLayes))) then
                                westN3coner_snd = westN3corner_snd + 1
                            elseif (i.le.xu .and. (i.gt.(xu-ghostLayes))) then
                                eastN3coner_snd = eastN3corner_snd + 1
                            endif

                            if (j.lt.yl) then
                                suthN3coner_rcv = suthN3corner_rcv + 1
                            elseif (j.gt.yu) then
                                nothN3coner_rcv = nothN3corner_rcv + 1
                            elseif (j.ge.yl .and. (j.lt.(yl+ghostLayes))) then
                                suthN3coner_snd = suthN3corner_snd + 1
                            elseif (j.le.yu .and. (j.gt.(yu-ghostLayes))) then
                                nothN3coner_snd = nothN3corner_snd + 1
                            endif

                            if (k.lt.zl) then
                                backN3coner_rcv = backN3corner_rcv + 1
                            elseif (k.gt.zu) then
                                fntN3corner_rcv = frntN3corner_rcv + 1
                            elseif (k.ge.zl .and. (k.lt.(zl+ghostLayes))) then
                                backN3coner_snd = backN3corner_snd + 1
                            elseif (k.le.zu .and. (k.gt.(zu-ghostLayes))) then
                                fntN3corner_snd = frntN3corner_snd + 1
                            endif
                        endif
                    endif
                enddo
            enddo
        enddo


        ! allocate fw only fo 3-fold corner, here nCorner counts only 3-fold
        ! coners
        allocate(fw(nConer, Nc8,1:8))

        ! corect the number of boundary processors' N3corner
        ! MPI_Waitall equire the snd and rcv buffer sizes to be the same
        ! even though the Y and Z diection is not periodical
        westGlbN3coner_rcv = 0
        eastGlbN3coner_rcv = 0
        suthGlbN3coner_rcv = 0
        nothGlbN3coner_rcv = 0
        backGlbN3coner_rcv = 0
        fntGlbN3corner_rcv = 0
        do k=zlg,zug
            do j=ylg,yug
                do i=1,ghostLayes
                    if(any(only3ConerLabels(:) == image(xl+i-1,j,k))) then
                        eastGlbN3coner_rcv = eastGlbN3corner_rcv+1
                    endif
                    if(any(only3ConerLabels(:) == image(xu-2+i,j,k))) then
                        westGlbN3coner_rcv = westGlbN3corner_rcv+1
                    endif
                enddo
            enddo
        enddo
        do k=zlg,zug
            do j=1,ghostLayes
                do i=xlg,xug
                    if(any(only3ConerLabels(:) == image(i,yl+j-1,k))) then
                        nothGlbN3coner_rcv = nothGlbN3corner_rcv+1
                    endif
                    if(any(only3ConerLabels(:) == image(i,yu-2+j,k))) then
                        suthGlbN3coner_rcv = suthGlbN3corner_rcv+1
                    endif
                enddo
            enddo
        enddo
        do k=1,ghostLayes
            do j=ylg,yug
                do i=xlg,xug
                    if(any(only3ConerLabels(:) == image(i,j,zl+k-1))) then
                        fntGlbN3corner_rcv = frntGlbN3corner_rcv+1
                    endif
                    if(any(only3ConerLabels(:) == image(i,j,zu-2+k))) then
                        backGlbN3coner_rcv = backGlbN3corner_rcv+1
                    endif
                enddo
            enddo
        enddo
        if(xl == xmin) westN3coner_rcv = westGlbN3corner_rcv
        if(xu == xmax) eastN3coner_rcv = eastGlbN3corner_rcv
        if(yl == ymin) suthN3coner_rcv = suthGlbN3corner_rcv
        if(yu == ymax) nothN3coner_rcv = nothGlbN3corner_rcv
        if(zl == zmin) backN3coner_rcv = backGlbN3corner_rcv
        if(zu == zmax) fntN3corner_rcv = frntGlbN3corner_rcv

        ! allocate maps fo 3-fold corners
        allocate(map3CoWsnd(westN3corner_snd))
        allocate(map3CoWrcv(westN3corner_rcv))
        allocate(map3CoEsnd(eastN3corner_snd))
        allocate(map3CoErcv(eastN3corner_rcv))
        allocate(map3CoSsnd(suthN3corner_snd))
        allocate(map3CoSrcv(suthN3corner_rcv))
        allocate(map3CoNsnd(nothN3corner_snd))
        allocate(map3CoNrcv(nothN3corner_rcv))
        allocate(map3CoBsnd(backN3corner_snd))
        allocate(map3CoBrcv(backN3corner_rcv))
        allocate(map3CoFsnd(frntN3corner_snd))
        allocate(map3CoFrcv(frntN3corner_rcv))

        ! constuct the maps for 3-fold corners
        eastN3coner_snd = 0
        eastN3coner_rcv = 0
        nothN3coner_snd = 0
        nothN3coner_rcv = 0
        fntN3corner_snd = 0
        fntN3corner_rcv = 0
        westN3coner_snd = 0
        westN3coner_rcv = 0
        suthN3coner_snd = 0
        suthN3coner_rcv = 0
        backN3coner_snd = 0
        backN3coner_rcv = 0
        nConer = 0
        do k=zlg, zug
            do j=ylg, yug
                do i=xlg, xug
                    localid = (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                    if(any(only3ConerLabels(:) == image(i,j,k))) then
                        nConer = nCorner + 1
                        if (i.lt.xl) then
                            westN3coner_rcv = westN3corner_rcv + 1
                            map3CoWrcv(westN3corner_rcv) = nCorner
                        elseif (i.gt.xu) then
                            eastN3coner_rcv = eastN3corner_rcv + 1
                            map3CoErcv(eastN3corner_rcv) = nCorner
                        elseif (i.ge.xl .and. (i.lt.(xl+ghostLayes))) then
                            westN3coner_snd = westN3corner_snd + 1
                            map3CoWsnd(westN3corner_snd) = nCorner
                        elseif (i.le.xu .and. (i.gt.(xu-ghostLayes))) then
                            eastN3coner_snd = eastN3corner_snd + 1
                            map3CoEsnd(eastN3corner_snd) = nCorner
                        endif

                        if (j.lt.yl) then
                            suthN3coner_rcv = suthN3corner_rcv + 1
                            map3CoSrcv(suthN3corner_rcv) = nCorner
                        elseif (j.gt.yu) then
                            nothN3coner_rcv = nothN3corner_rcv + 1
                            map3CoNrcv(nothN3corner_rcv) = nCorner
                        elseif (j.ge.yl .and. (j.lt.(yl+ghostLayes))) then
                            suthN3coner_snd = suthN3corner_snd + 1
                            map3CoSsnd(suthN3corner_snd) = nCorner
                        elseif (j.le.yu .and. (j.gt.(yu-ghostLayes))) then
                            nothN3coner_snd = nothN3corner_snd + 1
                            map3CoNsnd(nothN3corner_snd) = nCorner
                        endif

                        if (k.lt.zl) then
                            backN3coner_rcv = backN3corner_rcv + 1
                            map3CoBrcv(backN3corner_rcv) = nCorner
                        elseif (k.gt.zu) then
                            fntN3corner_rcv = frntN3corner_rcv + 1
                            map3CoFrcv(frntN3corner_rcv) = nCorner
                        elseif (k.ge.zl .and. (k.lt.(zl+ghostLayes))) then
                            backN3coner_snd = backN3corner_snd + 1
                            map3CoBsnd(backN3corner_snd) = nCorner
                        elseif (k.le.zu .and. (k.gt.(zu-ghostLayes))) then
                            fntN3corner_snd = frntN3corner_snd + 1
                            map3CoFsnd(frntN3corner_snd) = nCorner
                        endif
                    endif
                enddo
            enddo
        enddo
        !done constucting map



        ! allocate the aray for wall extrapolation coefficient
        allocate(extCoef(nWall,2))
        ! default 2nd oder
        extCoef(:, 1) =  2.d0
        extCoef(:, 2) = -1.d0
        do l=1, nWall
            localid = vecWall(l)
            k = zlg + localid/Nxytotal
            j = ylg + (localid - (k-zlg)*Nxytotal)/Nxtotal
            i = xlg + localid - (k-zlg)*Nxytotal - (j-ylg)*Nxtotal
            if(any(wallHasW(:) == image(i,j,k))) then
                if(i==xl+1 .and. image(xl,j,k) == fluid) then
                    extCoef(l,1) = 1.d0
                    extCoef(l,2) = 0.d0
                endif
            endif
            if(any(wallHasE(:) == image(i,j,k))) then
                if(i==xu-1 .and. image(xu,j,k) == fluid) then
                    extCoef(l,1) = 1.d0
                    extCoef(l,2) = 0.d0
                endif
            endif
            if(any(wallHasS(:) == image(i,j,k))) then
                if(j==yl+1 .and. image(i,yl,k) == fluid) then
                    extCoef(l,1) = 1.d0
                    extCoef(l,2) = 0.d0
                endif
            endif
            if(any(wallHasN(:) == image(i,j,k))) then
                if(j==yu-1 .and. image(i,yu,k) == fluid) then
                    extCoef(l,1) = 1.d0
                    extCoef(l,2) = 0.d0
                endif
            endif
            if(any(wallHasB(:) == image(i,j,k))) then
                if(k==zl+1 .and. image(i,j,zl) == fluid) then
                    extCoef(l,1) = 1.d0
                    extCoef(l,2) = 0.d0
                endif
            endif
            if(any(wallHasF(:) == image(i,j,k))) then
                if(k==zu-1 .and. image(i,j,zu) == fluid) then
                    extCoef(l,1) = 1.d0
                    extCoef(l,2) = 0.d0
                endif
            endif
        enddo

        if(wallExtOder == 2) then
            ! eset to 2nd order
            extCoef(:, 1) =  2.d0
            extCoef(:, 2) = -1.d0
        elseif (wallExtOder == 1) then
            ! eset to 1st order
            extCoef(:, 1) =  1.d0
            extCoef(:, 2) =  0.d0
        elseif (wallExtOder /= 3) then
            PRINT*, "Eror: wallExtOrder wroond shoud be [1|2|3]"
        endif

        ! count fluid points
        Nfluid = 0
        do k=zl,zu
            do j=yl,yu
                do i=xl,xu
                    if (image(i,j,k)==fluid) then
                        Nfluid=Nfluid+1
                    end if
                end do
            end do
        end do
        allocate(mapF(Nfluid))

        ! fill the map
        Nfluid = 0
        do k=zl,zu
            do j=yl,yu
                do i=xl,xu
                    if (image(i,j,k)==fluid) then
                        Nfluid=Nfluid+1
                        mapF(Nfluid) = (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                    end if
                end do
            end do
        end do

        !Diection 1
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xl == xmin) bxl = xl+1 !if most west block(pocessor)
        if(yl == ymin) byl = yl+1 !if most south block
        if(zl == zmin) bzl = zl+1!if most back  block
        Nstencil1=0
        do k=bzl,bzu
            do j=byl,byu
                do i=bxl,bxu
                    if (image(i,j,k)==fluid) then
                        Nstencil1=Nstencil1+1
                    end if
                end do
            end do
        end do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(di1(3, Nstencil1))
        icount=0
        do k=bzl,bzu
            do j=byl,byu
                do i=bxl,bxu
                    if (image(i,j,k)==fluid) then
                        icount=icount+1
                        di1(:,icount)=(/i,j,k/)
                    end if
                end do
            end do
        end do

        !Diection 2
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xu == xmax) bxu = xu-1 !if most west block(pocessor)
        if(yl == ymin) byl = yl+1 !if most south block
        if(zl == zmin) bzl = zl+1!if most back  block
        Nstencil2=0
        do k=bzl,bzu
            do j=byl,byu
                do i=bxu,bxl,-1
                    if (image(i,j,k)==fluid) then
                        Nstencil2=Nstencil2+1
                    end if
                end do
            end do
        end do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(di2(3, Nstencil2))
        icount=0
        do k=bzl,bzu
            do j=byl,byu
                do i=bxu,bxl,-1
                    if (image(i,j,k)==fluid) then
                        icount=icount+1
                        di2(:,icount)=(/i,j,k/)
                    end if
                end do
            end do
        end do

        !Diection 3
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xu == xmax) bxu = xu-1 !if most west block(pocessor)
        if(yu == ymax) byu = yu-1 !if most south block
        if(zl == zmin) bzl = zl+1!if most back  block
        Nstencil3=0
        do k=bzl,bzu
            do j=byu,byl,-1
                do i=bxu,bxl,-1
                    if (image(i,j,k)==fluid) then
                        Nstencil3=Nstencil3+1
                    end if
                end do
            end do
        end do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(di3(3, Nstencil3))
        icount=0
        do k=bzl,bzu
            do j=byu,byl,-1
                do i=bxu,bxl,-1
                    if (image(i,j,k)==fluid) then
                        icount=icount+1
                        di3(:,icount)=(/i,j,k/)
                    end if
                end do
            end do
        end do

        !Diection 4
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xl == xmin) bxl = xl+1 !if most west block(pocessor)
        if(yu == ymax) byu = yu-1 !if most south block
        if(zl == zmin) bzl = zl+1 !if most back  block
        Nstencil4=0
        do k=bzl,bzu
            do j=byu,byl,-1
                do i=bxl,bxu
                    if (image(i,j,k)==fluid) then
                        Nstencil4=Nstencil4+1
                    end if
                end do
            end do
        end do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(di4(3, Nstencil4))
        icount=0
        do k=bzl,bzu
            do j=byu,byl,-1
                do i=bxl,bxu
                    if (image(i,j,k)==fluid) then
                        icount=icount+1
                        di4(:,icount)=(/i,j,k/)
                    end if
                end do
            end do
        end do

        !Diection 5
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xl == xmin) bxl = xl+1 !if most west block(pocessor)
        if(yl == ymin) byl = yl+1 !if most south block
        if(zu == zmax) bzu = zu-1!if most back  block
        Nstencil5=0
        do k=bzu,bzl,-1
            do j=byl,byu
                do i=bxl,bxu
                    if (image(i,j,k)==fluid) then
                        Nstencil5=Nstencil5+1
                    end if
                end do
            end do
        end do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(di5(3, Nstencil5))
        icount=0
        do k=bzu,bzl,-1
            do j=byl,byu
                do i=bxl,bxu
                    if (image(i,j,k)==fluid) then
                        icount=icount+1
                        di5(:,icount)=(/i,j,k/)
                    end if
                end do
            end do
        end do

        !Diection 6
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xu == xmax) bxu = xu-1 !if most west block(pocessor)
        if(yl == ymin) byl = yl+1 !if most south block
        if(zu == zmax) bzu = zu-1!if most back  block
        Nstencil6=0
        do k=bzu,bzl,-1
            do j=byl,byu
                do i=bxu,bxl,-1
                    if (image(i,j,k)==fluid) then
                        Nstencil6=Nstencil6+1
                    end if
                end do
            end do
        end do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(di6(3, Nstencil6))
        icount=0
        do k=bzu,bzl,-1
            do j=byl,byu
                do i=bxu,bxl,-1
                    if (image(i,j,k)==fluid) then
                        icount=icount+1
                        di6(:,icount)=(/i,j,k/)
                    end if
                end do
            end do
        end do

        !Diection 7
        bxl = xl
        byl = yl
        bzl = zl
        bxu = xu
        byu = yu
        bzu = zu
        if(xu == xmax) bxu = xu-1 !if most west block(pocessor)
        if(yu == ymax) byu = yu-1 !if most south block
        if(zu == zmax) bzu = zu-1!if most back  block
        Nstencil7=0
        do k=bzu,bzl,-1
            do j=byu,byl,-1
                do i=bxu,bxl,-1
                    if (image(i,j,k)==fluid) then
                        Nstencil7=Nstencil7+1
                    end if
                end do
            end do
        end do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(di7(3, Nstencil7))
        icount=0
        do k=bzu,bzl,-1
            do j=byu,byl,-1
                do i=bxu,bxl,-1
                    if (image(i,j,k)==fluid) then
                        icount=icount+1
                        di7(:,icount)=(/i,j,k/)
                    end if
                end do
            end do
        end do

        !Diection 8
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
        do k=bzu,bzl,-1
            do j=byu,byl,-1
                do i=bxl,bxu
                    if (image(i,j,k)==fluid) then
                        Nstencil8=Nstencil8+1
                    end if
                end do
            end do
        end do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(di8(3, Nstencil8))
        icount=0
        do k=bzu,bzl,-1
            do j=byu,byl,-1
                do i=bxl,bxu
                    if (image(i,j,k)==fluid) then
                        icount=icount+1
                        di8(:,icount)=(/i,j,k/)
                    end if
                end do
            end do
        end do

        !constuct the deferencial coefficients
        allocate(coef1(Nstencil1,9),coef2(Nstencil2,9),&
            coef3(Nstencil3,9),coef4(Nstencil4,9))
        allocate(coef5(Nstencil5,9),coef6(Nstencil6,9),&
            coef7(Nstencil7,9),coef8(Nstencil8,9))
        
        do i=1,Nstencil1
            ii=di1(1,i) 
            jj=di1(2,i)
            kk=di1(3,i)
            !2nd oder of accuracy
            coef1(i,1)= 1.5d0/ds !x0n
            coef1(i,2)= 2.d0/ds  !x1n
            coef1(i,3)=-0.5d0/ds   !x2n
            coef1(i,4)= 1.5d0/ds !y0n
            coef1(i,5)= 2.d0/ds  !y1n
            coef1(i,6)=-0.5d0/ds   !y2n
            coef1(i,7)= 1.5d0/ds !z0n
            coef1(i,8)= 2.d0/ds  !z1n
            coef1(i,9)=-0.5d0/ds   !z2n
            
            !1st oder of accuracy in x
            if ((image(ii-2,jj,kk)==ghost)  .o.(image(ii-1,jj,kk)==wallE)  &
            .o.(image(ii-1,jj,kk)==wallES) .or.(image(ii-1,jj,kk)==wallEN) &
            .o.(image(ii-1,jj,kk)==wallEF) .or.(image(ii-1,jj,kk)==wallEB) &
            .o.(image(ii-1,jj,kk)==wallENF).or.(image(ii-1,jj,kk)==wallENB)&
            .o.(image(ii-1,jj,kk)==wallESF).or.(image(ii-1,jj,kk)==wallESB)) then
                coef1(i,1)=1.d0/ds  !x0n
                coef1(i,2)=1.d0/ds  !x1n
                coef1(i,3)=0.d0     !x2n
            end if
            !1st oder of accuracy in y
            if ((image(ii,jj-2,kk)==ghost)  .o.(image(ii,jj-1,kk)==wallN)  &
            .o.(image(ii,jj-1,kk)==wallEN) .or.(image(ii,jj-1,kk)==wallWN) &
            .o.(image(ii,jj-1,kk)==wallNF) .or.(image(ii,jj-1,kk)==wallNB) &
            .o.(image(ii,jj-1,kk)==wallENF).or.(image(ii,jj-1,kk)==wallENB)&
            .o.(image(ii,jj-1,kk)==wallWNF).or.(image(ii,jj-1,kk)==wallWNB)) then
                coef1(i,4)=1.d0/ds  !y0n
                coef1(i,5)=1.d0/ds  !y1n
                coef1(i,6)=0.d0     !y2n
            end if
            !1st oder of accuracy in z
            if ((image(ii,jj,kk-2)==ghost)  .o.(image(ii,jj,kk-1)==wallF)  &
            .o.(image(ii,jj,kk-1)==wallNF) .or.(image(ii,jj,kk-1)==wallSF) &
            .o.(image(ii,jj,kk-1)==wallEF) .or.(image(ii,jj,kk-1)==wallWF) &
            .o.(image(ii,jj,kk-1)==wallENF).or.(image(ii,jj,kk-1)==wallWNF)&
            .o.(image(ii,jj,kk-1)==wallESF).or.(image(ii,jj,kk-1)==wallWSF)) then
                coef1(i,7)=1.d0/ds  !z0n
                coef1(i,8)=1.d0/ds  !z1n
                coef1(i,9)=0.d0     !z2n
            end if
        end do ! end of ceof fo dir 1

        do i=1,Nstencil2
            ii=di2(1,i) 
            jj=di2(2,i)
            kk=di2(3,i)
            !2nd oder of accuracy
            coef2(i,1)=-1.5d0/ds !x0n
            coef2(i,2)=-2.d0/ds  !x1n
            coef2(i,3)=0.5d0/ds   !x2n
            coef2(i,4)=1.5d0/ds !y0n
            coef2(i,5)=2.d0/ds  !y1n
            coef2(i,6)=-0.5d0/ds   !y2n
            coef2(i,7)=1.5d0/ds !z0n
            coef2(i,8)=2.d0/ds  !z1n
            coef2(i,9)=-0.5d0/ds   !z2n
            !1st oder of accuracy in x
            if ((image(ii+2,jj,kk)==ghost)  .o.(image(ii+1,jj,kk)==wallW)  &
            .o.(image(ii+1,jj,kk)==wallWS) .or.(image(ii+1,jj,kk)==wallWN) &
            .o.(image(ii+1,jj,kk)==wallWF) .or.(image(ii+1,jj,kk)==wallWB) &
            .o.(image(ii+1,jj,kk)==wallWNF).or.(image(ii+1,jj,kk)==wallWNB)&
            .o.(image(ii+1,jj,kk)==wallWSF).or.(image(ii+1,jj,kk)==wallWSB)) then
                coef2(i,1)=-1.d0/ds  !x0n
                coef2(i,2)=-1.d0/ds  !x1n
                coef2(i,3)=0.d0     !x2n
            end if
            !1st oder of accuracy in y
            if ((image(ii,jj-2,kk)==ghost)  .o.(image(ii,jj-1,kk)==wallN)  &
            .o.(image(ii,jj-1,kk)==wallEN) .or.(image(ii,jj-1,kk)==wallWN) &
            .o.(image(ii,jj-1,kk)==wallNF) .or.(image(ii,jj-1,kk)==wallNB) &
            .o.(image(ii,jj-1,kk)==wallENF).or.(image(ii,jj-1,kk)==wallENB)&
            .o.(image(ii,jj-1,kk)==wallWNF).or.(image(ii,jj-1,kk)==wallWNB)) then
                coef2(i,4)=1.d0/ds  !y0n
                coef2(i,5)=1.d0/ds  !y1n
                coef2(i,6)=0.d0     !y2n
            end if
            !1st oder of accuracy in z
            if ((image(ii,jj,kk-2)==ghost)  .o.(image(ii,jj,kk-1)==wallF)  &
            .o.(image(ii,jj,kk-1)==wallNF) .or.(image(ii,jj,kk-1)==wallSF) &
            .o.(image(ii,jj,kk-1)==wallEF) .or.(image(ii,jj,kk-1)==wallWF) &
            .o.(image(ii,jj,kk-1)==wallENF).or.(image(ii,jj,kk-1)==wallWNF)&
            .o.(image(ii,jj,kk-1)==wallESF).or.(image(ii,jj,kk-1)==wallWSF)) then
                coef2(i,7)=1.d0/ds  !z0n
                coef2(i,8)=1.d0/ds  !z1n
                coef2(i,9)=0.d0     !z2n
            end if
        end do

        do i=1,Nstencil3
            ii=di3(1,i) 
            jj=di3(2,i)
            kk=di3(3,i)
            !2nd oder of accuracy
            coef3(i,1)=-1.5d0/ds !x0n
            coef3(i,2)=-2.d0/ds  !x1n
            coef3(i,3)=0.5d0/ds   !x2n
            coef3(i,4)=-1.5d0/ds !y0n
            coef3(i,5)=-2.d0/ds  !y1n
            coef3(i,6)=0.5d0/ds   !y2n
            coef3(i,7)=1.5d0/ds !z0n
            coef3(i,8)=2.d0/ds  !z1n
            coef3(i,9)=-0.5d0/ds   !z2n
            !1st oder of accuracy in x
            if ((image(ii+2,jj,kk)==ghost)  .o.(image(ii+1,jj,kk)==wallW)  &
            .o.(image(ii+1,jj,kk)==wallWS) .or.(image(ii+1,jj,kk)==wallWN) &
            .o.(image(ii+1,jj,kk)==wallWF) .or.(image(ii+1,jj,kk)==wallWB) &
            .o.(image(ii+1,jj,kk)==wallWNF).or.(image(ii+1,jj,kk)==wallWNB)&
            .o.(image(ii+1,jj,kk)==wallWSF).or.(image(ii+1,jj,kk)==wallWSB)) then
                coef3(i,1)=-1.d0/ds  !x0n
                coef3(i,2)=-1.d0/ds  !x1n
                coef3(i,3)=0.d0     !x2n
            end if
            !1st oder of accuracy in y
            if ((image(ii,jj+2,kk)==ghost)  .o.(image(ii,jj+1,kk)==wallS)  &
            .o.(image(ii,jj+1,kk)==wallES) .or.(image(ii,jj+1,kk)==wallWS) &
            .o.(image(ii,jj+1,kk)==wallSF) .or.(image(ii,jj+1,kk)==wallSB) &
            .o.(image(ii,jj+1,kk)==wallESF).or.(image(ii,jj+1,kk)==wallESB)&
            .o.(image(ii,jj+1,kk)==wallWSF).or.(image(ii,jj+1,kk)==wallWSB)) then
                coef3(i,4)=-1.d0/ds  !y0n
                coef3(i,5)=-1.d0/ds  !y1n
                coef3(i,6)=0.d0     !y2n
            end if
            !1st oder of accuracy in z
            if ((image(ii,jj,kk-2)==ghost)  .o.(image(ii,jj,kk-1)==wallF)  &
            .o.(image(ii,jj,kk-1)==wallNF) .or.(image(ii,jj,kk-1)==wallSF) &
            .o.(image(ii,jj,kk-1)==wallEF) .or.(image(ii,jj,kk-1)==wallWF) &
            .o.(image(ii,jj,kk-1)==wallENF).or.(image(ii,jj,kk-1)==wallWNF)&
            .o.(image(ii,jj,kk-1)==wallESF).or.(image(ii,jj,kk-1)==wallWSF)) then
                coef3(i,7)=1.d0/ds  !z0n
                coef3(i,8)=1.d0/ds  !z1n
                coef3(i,9)=0.d0     !z2n
            end if
        end do
        
        do i=1,Nstencil4
            ii=di4(1,i) 
            jj=di4(2,i)
            kk=di4(3,i)
            !2nd oder of accuracy
            coef4(i,1)=1.5d0/ds !x0n
            coef4(i,2)=2.d0/ds  !x1n
            coef4(i,3)=-0.5d0/ds   !x2n
            coef4(i,4)=-1.5d0/ds !y0n
            coef4(i,5)=-2.d0/ds  !y1n
            coef4(i,6)=0.5d0/ds   !y2n
            coef4(i,7)=1.5d0/ds !z0n
            coef4(i,8)=2.d0/ds  !z1n
            coef4(i,9)=-0.5d0/ds   !z2n
            !1st oder of accuracy in x
            if ((image(ii-2,jj,kk)==ghost)  .o.(image(ii-1,jj,kk)==wallE)   &
            .o.(image(ii-1,jj,kk)==wallES) .or.(image(ii-1,jj,kk)==wallEN)  &
            .o.(image(ii-1,jj,kk)==wallEF) .or.(image(ii-1,jj,kk)==wallEB)  &
            .o.(image(ii-1,jj,kk)==wallENF).or.(image(ii-1,jj,kk)==wallENB) &
            .o.(image(ii-1,jj,kk)==wallESF).or.(image(ii-1,jj,kk)==wallESB))  then
                coef4(i,1)=1.d0/ds  !x0n
                coef4(i,2)=1.d0/ds  !x1n
                coef4(i,3)=0.d0     !x2n
            end if
            !1st oder of accuracy in y
            if ((image(ii,jj+2,kk)==ghost)  .o.(image(ii,jj+1,kk)==wallS)  &
            .o.(image(ii,jj+1,kk)==wallES) .or.(image(ii,jj+1,kk)==wallWS) &
            .o.(image(ii,jj+1,kk)==wallSF) .or.(image(ii,jj+1,kk)==wallSB) &
            .o.(image(ii,jj+1,kk)==wallESF).or.(image(ii,jj+1,kk)==wallESB)&
            .o.(image(ii,jj+1,kk)==wallWSF).or.(image(ii,jj+1,kk)==wallWSB)) then
                coef4(i,4)=-1.d0/ds  !y0n
                coef4(i,5)=-1.d0/ds  !y1n
                coef4(i,6)=0.d0     !y2n
            end if
            !1st oder of accuracy in z
            if ((image(ii,jj,kk-2)==ghost)  .o.(image(ii,jj,kk-1)==wallF)  &
            .o.(image(ii,jj,kk-1)==wallNF) .or.(image(ii,jj,kk-1)==wallSF) &
            .o.(image(ii,jj,kk-1)==wallEF) .or.(image(ii,jj,kk-1)==wallWF) &
            .o.(image(ii,jj,kk-1)==wallENF).or.(image(ii,jj,kk-1)==wallWNF)&
            .o.(image(ii,jj,kk-1)==wallESF).or.(image(ii,jj,kk-1)==wallWSF)) then
                coef4(i,7)=1.d0/ds  !z0n
                coef4(i,8)=1.d0/ds  !z1n
                coef4(i,9)=0.d0     !z2n
            end if
        end do
        
        do i=1,Nstencil5
            ii=di5(1,i) 
            jj=di5(2,i)
            kk=di5(3,i)
            !2nd oder of accuracy
            coef5(i,1)=1.5d0/ds !x0n
            coef5(i,2)=2.d0/ds  !x1n
            coef5(i,3)=-0.5d0/ds   !x2n
            coef5(i,4)=1.5d0/ds !y0n
            coef5(i,5)=2.d0/ds  !y1n
            coef5(i,6)=-0.5d0/ds   !y2n
            coef5(i,7)=-1.5d0/ds !z0n
            coef5(i,8)=-2.d0/ds  !z1n
            coef5(i,9)=0.5d0/ds   !z2n
            !1st oder of accuracy in x
            if ((image(ii-2,jj,kk)==ghost)  .o.(image(ii-1,jj,kk)==wallE)  &
            .o.(image(ii-1,jj,kk)==wallES) .or.(image(ii-1,jj,kk)==wallEN) &
            .o.(image(ii-1,jj,kk)==wallEF) .or.(image(ii-1,jj,kk)==wallEB) &
            .o.(image(ii-1,jj,kk)==wallENF).or.(image(ii-1,jj,kk)==wallENB)&
            .o.(image(ii-1,jj,kk)==wallESF).or.(image(ii-1,jj,kk)==wallESB)) then
                coef5(i,1)=1.d0/ds  !x0n
                coef5(i,2)=1.d0/ds  !x1n
                coef5(i,3)=0.d0     !x2n
            end if
            !1st oder of accuracy in y
            if ((image(ii,jj-2,kk)==ghost)  .o.(image(ii,jj-1,kk)==wallN)  &
            .o.(image(ii,jj-1,kk)==wallEN) .or.(image(ii,jj-1,kk)==wallWN) &
            .o.(image(ii,jj-1,kk)==wallNF) .or.(image(ii,jj-1,kk)==wallNB) &
            .o.(image(ii,jj-1,kk)==wallENF).or.(image(ii,jj-1,kk)==wallENB)&
            .o.(image(ii,jj-1,kk)==wallWNF).or.(image(ii,jj-1,kk)==wallWNB)) then
                coef5(i,4)=1.d0/ds  !y0n
                coef5(i,5)=1.d0/ds  !y1n
                coef5(i,6)=0.d0     !y2n
            end if
            !1st oder of accuracy in z
            if ((image(ii,jj,kk+2)==ghost)  .o.(image(ii,jj,kk+1)==wallB)  &
            .o.(image(ii,jj,kk+1)==wallNB) .or.(image(ii,jj,kk+1)==wallSB) &
            .o.(image(ii,jj,kk+1)==wallEB) .or.(image(ii,jj,kk+1)==wallWB) &
            .o.(image(ii,jj,kk+1)==wallENB).or.(image(ii,jj,kk+1)==wallWNB)&
            .o.(image(ii,jj,kk+1)==wallESB).or.(image(ii,jj,kk+1)==wallWSB)) then
                coef5(i,7)=-1.d0/ds  !z0n
                coef5(i,8)=-1.d0/ds  !z1n
                coef5(i,9)=0.d0     !z2n
            end if
        end do

        do i=1,Nstencil6
            ii=di6(1,i) 
            jj=di6(2,i)
            kk=di6(3,i)
            !2nd oder of accuracy
            coef6(i,1)=-1.5d0/ds !x0n
            coef6(i,2)=-2.d0/ds  !x1n
            coef6(i,3)=0.5d0/ds   !x2n
            coef6(i,4)=1.5d0/ds !y0n
            coef6(i,5)=2.d0/ds  !y1n
            coef6(i,6)=-0.5d0/ds   !y2n
            coef6(i,7)=-1.5d0/ds !z0n
            coef6(i,8)=-2.d0/ds  !z1n
            coef6(i,9)=0.5d0/ds   !z2n
            !1st oder of accuracy in x
            if ((image(ii+2,jj,kk)==ghost)  .o.(image(ii+1,jj,kk)==wallW)  &
            .o.(image(ii+1,jj,kk)==wallWS) .or.(image(ii+1,jj,kk)==wallWN) &
            .o.(image(ii+1,jj,kk)==wallWF) .or.(image(ii+1,jj,kk)==wallWB) &
            .o.(image(ii+1,jj,kk)==wallWNF).or.(image(ii+1,jj,kk)==wallWNB)&
            .o.(image(ii+1,jj,kk)==wallWSF).or.(image(ii+1,jj,kk)==wallWSB)) then
                coef6(i,1)=-1.d0/ds  !x0n
                coef6(i,2)=-1.d0/ds  !x1n
                coef6(i,3)=0.d0     !x2n
            end if
            !1st oder of accuracy in y
            if ((image(ii,jj-2,kk)==ghost)  .o.(image(ii,jj-1,kk)==wallN)  &
            .o.(image(ii,jj-1,kk)==wallEN) .or.(image(ii,jj-1,kk)==wallWN) &
            .o.(image(ii,jj-1,kk)==wallNF) .or.(image(ii,jj-1,kk)==wallNB) &
            .o.(image(ii,jj-1,kk)==wallENF).or.(image(ii,jj-1,kk)==wallENB)&
            .o.(image(ii,jj-1,kk)==wallWNF).or.(image(ii,jj-1,kk)==wallWNB))  then
                coef6(i,4)=1.d0/ds  !y0n
                coef6(i,5)=1.d0/ds  !y1n
                coef6(i,6)=0.d0     !y2n
            end if
            !1st oder of accuracy in z
            if ((image(ii,jj,kk+2)==ghost)  .o.(image(ii,jj,kk+1)==wallB)  &
            .o.(image(ii,jj,kk+1)==wallNB) .or.(image(ii,jj,kk+1)==wallSB) &
            .o.(image(ii,jj,kk+1)==wallEB) .or.(image(ii,jj,kk+1)==wallWB) &
            .o.(image(ii,jj,kk+1)==wallENB).or.(image(ii,jj,kk+1)==wallWNB)&
            .o.(image(ii,jj,kk+1)==wallESB).or.(image(ii,jj,kk+1)==wallWSB)) then
                coef6(i,7)=-1.d0/ds  !z0n
                coef6(i,8)=-1.d0/ds  !z1n
                coef6(i,9)=0.d0     !z2n
            end if
        end do

        do i=1,Nstencil7
            ii=di7(1,i) 
            jj=di7(2,i)
            kk=di7(3,i)
            !2nd oder of accuracy
            coef7(i,1)=-1.5d0/ds !x0n
            coef7(i,2)=-2.d0/ds  !x1n
            coef7(i,3)=0.5d0/ds   !x2n
            coef7(i,4)=-1.5d0/ds !y0n
            coef7(i,5)=-2.d0/ds  !y1n
            coef7(i,6)=0.5d0/ds   !y2
            coef7(i,7)=-1.5d0/ds !z0n
            coef7(i,8)=-2.d0/ds  !z1n
            coef7(i,9)=0.5d0/ds   !z2n
            !1st oder of accuracy in x
            if ((image(ii+2,jj,kk)==ghost)  .o.(image(ii+1,jj,kk)==wallW)  &
            .o.(image(ii+1,jj,kk)==wallWS) .or.(image(ii+1,jj,kk)==wallWN) &
            .o.(image(ii+1,jj,kk)==wallWF) .or.(image(ii+1,jj,kk)==wallWB) &
            .o.(image(ii+1,jj,kk)==wallWNF).or.(image(ii+1,jj,kk)==wallWNB)&
            .o.(image(ii+1,jj,kk)==wallWSF).or.(image(ii+1,jj,kk)==wallWSB)) then
                coef7(i,1)=-1.d0/ds  !x0n
                coef7(i,2)=-1.d0/ds  !x1n
                coef7(i,3)=0.d0     !x2n
            end if
            !1st oder of accuracy in y
            if ((image(ii,jj+2,kk)==ghost)  .o.(image(ii,jj+1,kk)==wallS)  &
            .o.(image(ii,jj+1,kk)==wallES) .or.(image(ii,jj+1,kk)==wallWS) &
            .o.(image(ii,jj+1,kk)==wallSF) .or.(image(ii,jj+1,kk)==wallSB) &
            .o.(image(ii,jj+1,kk)==wallESF).or.(image(ii,jj+1,kk)==wallESB)&
            .o.(image(ii,jj+1,kk)==wallWSF).or.(image(ii,jj+1,kk)==wallWSB)) then
                coef7(i,4)=-1.d0/ds  !y0n
                coef7(i,5)=-1.d0/ds  !y1n
                coef7(i,6)=0.d0     !y2n
            end if
            !1st oder of accuracy in z
            if ((image(ii,jj,kk+2)==ghost)  .o.(image(ii,jj,kk+1)==wallB)  &
            .o.(image(ii,jj,kk+1)==wallNB) .or.(image(ii,jj,kk+1)==wallSB) &
            .o.(image(ii,jj,kk+1)==wallEB) .or.(image(ii,jj,kk+1)==wallWB) &
            .o.(image(ii,jj,kk+1)==wallENB).or.(image(ii,jj,kk+1)==wallWNB)&
            .o.(image(ii,jj,kk+1)==wallESB).or.(image(ii,jj,kk+1)==wallWSB)) then
                coef7(i,7)=-1.d0/ds  !z0n
                coef7(i,8)=-1.d0/ds  !z1n
                coef7(i,9)=0.d0     !z2n
            end if
        end do


        do i=1,Nstencil8
            ii=di8(1,i) 
            jj=di8(2,i)
            kk=di8(3,i)
            !2nd oder of accuracy
            coef8(i,1)=1.5d0/ds !x0n
            coef8(i,2)=2.d0/ds  !x1n
            coef8(i,3)=-0.5d0/ds   !x2n
            coef8(i,4)=-1.5d0/ds !y0n
            coef8(i,5)=-2.d0/ds  !y1n
            coef8(i,6)=0.5d0/ds   !y2n
            coef8(i,7)=-1.5d0/ds !z0n
            coef8(i,8)=-2.d0/ds  !z1n
            coef8(i,9)=0.5d0/ds   !z2n
            !1st oder of accuracy in x
            if ((image(ii-2,jj,kk)==ghost)  .o.(image(ii-1,jj,kk)==wallE)  &
            .o.(image(ii-1,jj,kk)==wallES) .or.(image(ii-1,jj,kk)==wallEN) &
            .o.(image(ii-1,jj,kk)==wallEF) .or.(image(ii-1,jj,kk)==wallEB) &
            .o.(image(ii-1,jj,kk)==wallENF).or.(image(ii-1,jj,kk)==wallENB)&
            .o.(image(ii-1,jj,kk)==wallESF).or.(image(ii-1,jj,kk)==wallESB)) then
                coef8(i,1)=1.d0/ds  !x0n
                coef8(i,2)=1.d0/ds  !x1n
                coef8(i,3)=0.d0     !x2n
            end if
            !1st oder of accuracy in y
            if ((image(ii,jj+2,kk)==ghost)  .o.(image(ii,jj+1,kk)==wallS)  &
            .o.(image(ii,jj+1,kk)==wallES) .or.(image(ii,jj+1,kk)==wallWS) &
            .o.(image(ii,jj+1,kk)==wallSF) .or.(image(ii,jj+1,kk)==wallSB) &
            .o.(image(ii,jj+1,kk)==wallESF).or.(image(ii,jj+1,kk)==wallESB)&
            .o.(image(ii,jj+1,kk)==wallWSF).or.(image(ii,jj+1,kk)==wallWSB)) then
                coef8(i,4)=-1.d0/ds  !y0n
                coef8(i,5)=-1.d0/ds  !y1n
                coef8(i,6)=0.d0     !y2n
            end if
            !1st oder of accuracy in z
            if ((image(ii,jj,kk+2)==ghost  ).o.(image(ii,jj,kk+1)==wallB)  &
            .o.(image(ii,jj,kk+1)==wallNB ).or.(image(ii,jj,kk+1)==wallSB) &
            .o.(image(ii,jj,kk+1)==wallEB ).or.(image(ii,jj,kk+1)==wallWB) &
            .o.(image(ii,jj,kk+1)==wallENB).or.(image(ii,jj,kk+1)==wallWNB)&
            .o.(image(ii,jj,kk+1)==wallESB).or.(image(ii,jj,kk+1)==wallWSB)) then
                coef8(i,7)=-1.d0/ds  !z0n
                coef8(i,8)=-1.d0/ds  !z1n
                coef8(i,9)=0.d0     !z2n
            end if
        end do

        fw=0.d0
    end suboutine setupPhysicalGrid
end module physicalGid
