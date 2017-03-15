!=======================================================================
!> @brief Physical space configurations
!=======================================================================
module physicalGrid
IMPLICIT NONE
SAVE

! Domain size

integer :: xl, xu, yl, yu, zl, zu
integer :: xlg, xug, ylg, yug, zlg, zug
integer, parameter :: ghostLayers = 2
integer, parameter :: fluidlayer=4
integer, parameter :: translate=0

!Obstacle geometry
integer, parameter :: obstR = floor((Nx-1)*((1.d0-porosity)*3.d0/4.d0/PI)**(1.d0/3.d0))
integer, parameter :: obstX = ghostLayer+Ny-2   !Check
integer, parameter :: obstY = ghostLayer-1      !Check
integer, parameter :: obstZ = ghostLayer-1  

! NX and NY is the global grid size
!integer, parameter :: Nx = 533, Ny = 428 !Brea stone
integer, parameter :: Ny=400
integer, parameter :: Nx=Ny+2*fluidlayer
integer, parameter :: Nz=Ny

integer, parameter :: xmin = 1
integer, parameter :: xmax = Nx
integer, parameter :: ymin = 1
integer, parameter :: ymax = Ny
integer, parameter :: zmin = 1
integer, parameter :: zmax = Nz

! need to be determined from mpi cood
integer :: Nxtotal, Nytotal, Nztotal, Nzytotal, Nxsub, Nysub, Nzsub, Ntotal

double precision, parameter :: ds =  1.0d0/(Nx-1)
integer, parameter :: column=fluidlayer/2 ! layer to extract flow rate

! raw and extended flag array
integer, dimension (:,:,:), allocatable :: array3D, array3Dg ! globle image
integer, dimension (:,:,:), allocatable :: image ! local image 

! grid point flags and wall type id
integer, parameter:: fluid = 0, solid = 1, ghost=4
integer, parameter:: wallE = 20, wallW = 21, wallN = 22, wallS = 23, wallF = 24, wallB = 25 !(1-side wall )
integer, parameter:: wallEN = 30, wallWN = 31, wallES = 32, wallWS = 33 !(2-side wall)
integer, parameter:: wallNF = 40, wallNB = 41, wallSF = 42, wallSB = 43 !(2-side wall)
integer, parameter:: wallEF = 50, wallEB = 51, wallWF = 52, wallWB = 53 !(2-side wall)
integer, parameter:: wallENF = 60, wallWNF = 61, wallWSF = 62, wallESF = 63 !(3-side wall)
integer, parameter:: wallENB = 70, wallWNB = 71, wallWSB = 72, wallESB = 73 !(3-side wall)

integer :: Nstencil1, Nstencil2, Nstencil3, Nstencil4 !(group 1-4)
integer :: Nstencil5, Nstencil6, Nstencil7, Nstencil8 !(group 5-8)

double precision :: real_porosity

! number of wall points
integer :: nWall 

integer, DIMENSION(:), ALLOCATABLE :: vecWall
integer, DIMENSION(:), ALLOCATABLE :: dir1, dir2, dir3, dir4, dir5, dir6, dir7, dir8! dir1(i) is the ith 
double precision, DIMENSION(:,:), ALLOCATABLE :: coef1, coef2, coef3, coef4, coef5, coef6, coef7, coef8
double precision, DIMENSION(:,:), ALLOCATABLE :: f1wZ,f2wZ,f3wZ,f4wZ,f5wZ,f6wZ,f7wZ,f8wZ


contains 
    ! should be called after calling MPIParams::setupVirtualProcessGrid
    ! such that xl, xlg, xu, xug and yl, ylg, yu, yug has been set already
    subroutine setupPhysicalGrid
        implicit none

        ! local vars
        integer :: localid, i, j, k, icount, countVoidP, NneighborFluid
        integer :: bxl, bxu, byl, byu, bzl, bzu ! bound when counting fluid point for sweeping
        integer :: ii, jj, kk
        character(kind=1) :: ctemp

        ! set the extend and sizes
        Nxtotal = xug - xlg + 1
        Nytotal = yug - ylg + 1
        Nytotal = zug - zlg + 1
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

        allocate(image(xlg:xug, ylg:yug, zlg:zug)) ! this is local one


        ! set global image array
        array3D  = fluid ! set all to fluid first (fluidLayers are fluid)
        array3Dg = fluid ! set to ghost, then set inner ridge

        ! set ghost
        do k=zmin-ghostLayers,zmax+ghostLayers
            do j=ymin-ghostLayers,ymax+ghostLayers
                do i=xmin-ghostLayers,xmax+ghostLayers
                    if(k<zmin .or. k>zmax .or. y<ymin .or. y>ymax .or. x<xmin .or. x>xmax) then
                        array3Dg(i,j,k) = ghost
                    endif
                enddo
            enddo
        enddo

        !read digital image, NOTE that all raw image dimension are NY^3. 
        Open(200,file='shale400^3.raw',status='OLD',form='unformatted',ACCESS="STREAM")
        Do k=1,Nz
            Do j=1,Ny
                Do i=fluidlayer+1,Ny+fluidlayer  
                    read(200) ctemp       
                    !l=(i-translate+ghostLayer)+(j-translate+ghostLayer-1)*Nxtotal+(k-translate+ghostLayer-1)*Nxytotal
                    if ((k>=1+translate).AND.(k<=Nz+translate).AND.(j>=1+translate).AND.(j<=Ny+translate).AND.&
                        & (i>=fluidlayer+1+translate).AND.(i<=Nx-fluidlayer+translate)) then
                        if (ichar(ctemp)==0) then
                            array3D(i,j,k) = fluid
                            array3Dg(i,j,k) = fluid
                        else if (ichar(ctemp)==1) then
                            array3D(i,j,k) = solid
                            array3Dg(i,j,k) = solid
                        else 
                            write(*,*) "invalid image data file"
                            stop 0
                        endif
                    endif
                Enddo  
            Enddo
        Enddo
        Close(200)

        ! drill the small pores (change array3Dg)
        do k=1,Nz
            do j=1,Ny
                do i=2,Nx-1
                    if(array3Dg(i,j,k) == fluid) then
                        if (array3Dg(i+1,j,k) /=fluid .and. array3Dg(i-1,j,k) /= fluid .or.
                            array3Dg(i,j+1,k) /=fluid .and. array3Dg(i,j-1,k) /= fluid .or.
                            array3Dg(i,j,k+1) /=fluid .and. array3Dg(i,j,k-1)) then
                            if (j==Ny) then
                                if (k==Nz) then
                                    array3Dg(i:i+1,j-1:j,k-1:k) = fluid
                                else
                                    array3Dg(i:i+1,j-1:j,k:k+1) = fluid
                                endif
                            else if (k==Nz) then
                                array3Dg(i:i+1,j:j+1,k-1,k) = fluid
                            else
                                array3Dg(i:i+1,j:j+1,k:k+1) = fluid
                            endif
                        endif
                    endif !fluid
                enddo !i
            enddo !j
        enddo !k

        ! set array3Dg
        array3Dg = ghost
        do k = zmin, zmax
            do j = ymin, ymax
                do i = xmin, xmax
                    array2Dg(i,j) = array2D(i,j)
                enddo
            enddo
        enddo


        ! set local image


        ! set local image
        image = fluid
        bxl = xlg
        byl = ylg
        bxu = xug
        byu = yug
        if(xl == xmin) bxl = xl !if most west block(processor)
        if(xu == xmax) bxu = xu !if most east block(processor)
        if(yl == ymin) byl = yl !if most south block
        if(yu == ymax) byu = yu !if most north block
        do j = byl, byu
            do i = bxl, bxu
               localid = (j-ylg)*Nxtotal + i-xlg+1
               image(localid) = array2D(i,j)
           enddo
        enddo

        !Assign the outerboarder layers (called ghost point in serial program)
        If (ghostLayers>0) then
            Do j=ylg,yug
                Do i=xlg,xug
                    ! test if boundary processor, here the ghost flag doesn't means 
                    ! the communication boundary, but only means the true global 
                    ! outter boarder.
                    If ((i<xmin).OR.(i>xmax).OR.(j<ymin).OR.(j>ymax)) then 
                        localid = (j-ylg)*Nxtotal + i-xlg+1
                        image(localid) = ghost
                    End if
                Enddo
            Enddo
        End if

        bxl = xlg
        bxu = xug
        byl = ylg
        byu = yug
        if(xl == xmin) bxl = xl !if most west block(processor)
        if(xu == xmax) bxu = xu !if most east block(processor)
        if(yl == ymin) byl = yl !if most south block
        if(yu == ymax) byu = yu !if most north block

        ! set wall points type based on sournding point type(f/s)
        nWall=0 ! count the wall points
        !print*, bxl, bxu, byl, byu
        Do j=byl,byu
            Do i=bxl,bxu
                !ii = i-xl+1
                !jj = j-yl+1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (array2D(i,j)==solid) then
                    NneighborFluid=1 !(1-2-3-4 in D2Q9 corresponding to 2-3-5-7)
                    ! found bug here, array index out of bound of array2D
                    !If (image(localid+1)==fluid)  NneighborFluid=NneighborFluid*2 
                    !If (image(localid+Nxtotal)==fluid)  NneighborFluid=NneighborFluid*3
                    !If (image(localid-1)==fluid)  NneighborFluid=NneighborFluid*5 
                    !If (image(localid-Nxtotal)==fluid)  NneighborFluid=NneighborFluid*7
                    If (array2Dg(i+1, j)==fluid)  NneighborFluid=NneighborFluid*2   !Check neighbor on the East(2)
                    If (array2Dg(i, j+1)==fluid)  NneighborFluid=NneighborFluid*3   !Check neighbor on the North(3)
                    If (array2Dg(i-1, j)==fluid)  NneighborFluid=NneighborFluid*5   !Check neighbor on the West(5)
                    If (array2Dg(i, j-1)==fluid)  NneighborFluid=NneighborFluid*7   !Check neighbor on the South(7)

                    SELECT Case (NneighborFluid)
                        CASE (2)
                            image(localid)=WallXp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (3)
                            image(localid)=WallYp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (5)
                            image(localid)=WallXn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (7)
                            image(localid)=WallYn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (6)
                            image(localid)=WallXpYp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (15)
                            image(localid)=WallXnYp
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (35)
                            image(localid)=WallXnYn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                        CASE (14)
                            image(localid)=WallXpYn
                            if(j .ge. yl .and. j .le. yu .and. i .ge. xl .and. i .le. xu) nWall=nWall+1
                    END SELECT
                Endif
            Enddo
        Enddo

        PRINT*, "nWall = ", nWall

        ! Create wall-type vectors
        !vecWall(walli) mark the global id of the walli'th wall in the image 
        ALLOCATE(vecWall(nWall))
        nWall=0
        Do j=yl, yu
            Do i=xl, xu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                SELECT Case (image(localid))
                    CASE (WallXp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallYp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallYn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXpYp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXnYp)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXnYn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                    CASE (WallXpYn)
                        nWall=nWall+1
                        vecWall(nWall)=localid
                END SELECT
            Enddo
        Enddo

        !BUG FOUND
        !Direction 1
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if(xl == xmin) bxl = xl + 1 !if most west block(processor)
        if(yl == ymin) byl = yl + 1 !if most south block
        !if(xu == xmax) bxu = xu     !if most east block(processor)
        !if(yu == ymax) byu = yu     !if most north block   
        !count fluid points when sweeping from 1st direction
        Nstencil1=0
        Do j=byl,byu
            Do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    Nstencil1=Nstencil1+1
                End if
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir1(Nstencil1))
        icount=0
        Do j=byl,byu
            Do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    icount=icount+1
                    dir1(icount)=localid
                End if
            End do
        End do

        !Direction 2
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if(xu == xmax) bxu = xu - 1 !if most east block
        if(yl == ymin) byl = yl + 1 !if most south block
        !if(xl == xmin) bxl = xl     !if most west block
        !if(yu == ymax) byu = yu     !if most north block
        !count fluid points when sweeping from 2nd direction
        Nstencil2=0
        Do j=byl,byu
            Do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    Nstencil2=Nstencil2+1
                End if
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir2(Nstencil2))
        icount=0
        Do j=byl,byu
            Do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    icount=icount+1
                    dir2(icount)=localid
                End if
            End do
        End do

        !Direction 3
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if(xu == xmax) bxu = xu - 1 !if most east block
        if(yu == ymax) byu = yu - 1 !if most north block
        !if(xl == xmin) bxl = xl     !if most west block
        !if(yl == ymin) byl = yl     !if most south block
        !count fluid points when sweeping from 3rd direction
        Nstencil3=0
        Do j=byu,byl,-1
            Do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    Nstencil3=Nstencil3+1
                End if
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir3(Nstencil3))
        icount=0
        Do j=byu,byl,-1
            Do i=bxu,bxl,-1
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    icount=icount+1
                    dir3(icount)=localid
                End if
            End do
        End do

        !Direction 4
        !bound
        bxl = xl
        bxu = xu
        byl = yl
        byu = yu
        if(xl == xmin) bxl = xl + 1 !if most west block
        if(yu == ymax) byu = yu - 1 !if most north block
        !if(xu == xmax) bxu = xu     !if most east block
        !if(yl == ymin) byl = yl     !if most south block
        !count fluid points when sweeping from 4th direction
        Nstencil4=0
        Do j=byu,byl,-1
            Do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    Nstencil4=Nstencil4+1
                End if
            End do
        End do
        !set the icount'th fluid node's localid in the 2D patch
        allocate(dir4(Nstencil4))
        icount=0
        Do j=byu,byl,-1
            Do i=bxl,bxu
                localid = (j-ylg)*Nxtotal + i-xlg+1
                If (image(localid)==fluid) then
                    icount=icount+1
                    dir4(icount)=localid
                End if
            End do
        End do

        !construct the deferencial coefficients
        ALLOCATE(coefI(Nstencil1,6),coefII(Nstencil2,6),coefIII(Nstencil3,6),coefIV(Nstencil4,6))
        Do i=1,Nstencil1
            localid=dir1(i) 
            !2nd order of accuracy
            coefI(i,1)=1.5d0/ds !x0n
            coefI(i,2)=2.d0/ds  !x1n
            coefI(i,3)=-0.5d0/ds   !x2n
            coefI(i,4)=1.5d0/ds !y0n
            coefI(i,5)=2.d0/ds  !y1n
            coefI(i,6)=-0.5d0/ds   !y2n
            if ((image(localid-2)==ghost).OR.(image(localid-1)==WallXp).OR.(image(localid-1)==WallXpYn) &
            & .OR.(image(localid-1)==WallXpYp)) then !1st order of accuracy in x
                coefI(i,1)=1.d0/ds  !x0n
                coefI(i,2)=1.d0/ds  !x1n
                coefI(i,3)=0.d0     !x2n
            end if
            if ((image(localid-2*Nxtotal)==ghost).OR.(image(localid-Nxtotal)==WallYp) &
            & .OR.(image(localid-Nxtotal)==WallXpYp).OR.(image(localid-Nxtotal)==WallXnYp)) then   !1st order of accuracy in y
                coefI(i,4)=1.d0/ds  !y0n
                coefI(i,5)=1.d0/ds  !y1n
                coefI(i,6)=0.d0     !y2n
            end if
        End do

        Do i=1,Nstencil2
            localid=dir2(i)
            !2nd order of accuracy
            coefII(i,1)=-1.5d0/ds !x0n
            coefII(i,2)=-2.d0/ds  !x1n
            coefII(i,3)=0.5d0/ds   !x2n
            coefII(i,4)=1.5d0/ds !y0n
            coefII(i,5)=2.d0/ds  !y1n
            coefII(i,6)=-0.5d0/ds   !y2n
            if ((image(localid+2)==ghost).OR.(image(localid+1)==WallXn).OR.(image(localid+1)==WallXnYp) &
            & .OR.(image(localid+1)==WallXnYn)) then !1st order of accuracy in x
                coefII(i,1)=-1.d0/ds  !x0n
                coefII(i,2)=-1.d0/ds  !x1n
                coefII(i,3)=0.d0     !x2n
            end if
            if ((image(localid-2*Nxtotal)==ghost).OR.(image(localid-Nxtotal)==WallYp) &
            & .OR.(image(localid-Nxtotal)==WallXpYp).OR.(image(localid-Nxtotal)==WallXnYp)) then   !1st order of accuracy in y
                coefII(i,4)=1.d0/ds  !y0n
                coefII(i,5)=1.d0/ds  !y1n
                coefII(i,6)=0.d0     !y2n
            end if
        End do

        Do i=1,Nstencil3
            localid=dir3(i)
            !2nd order of accuracy
            coefIII(i,1)=-1.5d0/ds !x0n
            coefIII(i,2)=-2.d0/ds  !x1n
            coefIII(i,3)=0.5d0/ds   !x2n
            coefIII(i,4)=-1.5d0/ds !y0n
            coefIII(i,5)=-2.d0/ds  !y1n
            coefIII(i,6)=0.5d0/ds   !y2n
            if ((image(localid+2)==ghost).OR.(image(localid+1)==WallXn).OR.(image(localid+1)==WallXnYp) &
            & .OR.(image(localid+1)==WallXnYn)) then !1st order of accuracy in x
                coefIII(i,1)=-1.d0/ds  !x0n
                coefIII(i,2)=-1.d0/ds  !x1n
                coefIII(i,3)=0.d0     !x2n
            end if
            if ((image(localid+2*Nxtotal)==ghost).OR.(image(localid+Nxtotal)==WallYn) &
            & .OR.(image(localid+Nxtotal)==WallXnYn).OR.(image(localid+Nxtotal)==WallXpYn)) then   !1st order of accuracy in y
                coefIII(i,4)=-1.d0/ds  !y0n
                coefIII(i,5)=-1.d0/ds  !y1n
                coefIII(i,6)=0.d0     !y2n
            end if
        End do

        Do i=1,Nstencil4
            localid=dir4(i)
            !2nd order of accuracy
            coefIV(i,1)=1.5d0/ds !x0n
            coefIV(i,2)=2.d0/ds  !x1n
            coefIV(i,3)=-0.5d0/ds   !x2n
            coefIV(i,4)=-1.5d0/ds !y0n
            coefIV(i,5)=-2.d0/ds  !y1n
            coefIV(i,6)=0.5d0/ds   !y2n
            if ((image(localid-2)==ghost).OR.(image(localid-1)==WallXp).OR.(image(localid-1)==WallXpYn) &
            & .OR.(image(localid-1)==WallXpYp)) then !1st order of accuracy in x
                coefIV(i,1)=1.d0/ds  !x0n
                coefIV(i,2)=1.d0/ds  !x1n
                coefIV(i,3)=0.d0     !x2n
            end if
            if ((image(localid+2*Nxtotal)==ghost).OR.(image(localid+Nxtotal)==WallYn) &
            & .OR.(image(localid+Nxtotal)==WallXnYn).OR.(image(localid+Nxtotal)==WallXpYn)) then   !1st order of accuracy in y
                coefIV(i,4)=-1.d0/ds  !y0n
                coefIV(i,5)=-1.d0/ds  !y1n
                coefIV(i,6)=0.d0     !y2n
            end if
        End do
    end subroutine setupPhysicalGrid
end module physicalGrid
