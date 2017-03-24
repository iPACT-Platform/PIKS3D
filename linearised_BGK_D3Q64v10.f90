!
!  This program calculates the pressure driven gas flow throught a
!  three-dimensional channel with obstacle using the linearised BGK kinetic equation
!  fix-point method
!
!
!       ---------------------------------
!       ---------------------------------
!       ||i                           o||                   y
!       ||n                           u||                   /\
!       ||-                           t||                   |
!       ||l                           l||                   |
!       ||e                           e||                   |
!       ||t                           t||                   |
!       ---------------------------------                   |----------->x
!       ---------------------------------                  /
!                                                         /
!                                                        /
!                                                     z /
!
!
!
!                               |
!                      WNF/B    |       ENF/B
!                               |
!               ----------------------------------  8 molecular velocity groups
!                               |
!                      WSF/B    |       ESF/B
!                               |
! Symmetry on lateral wall
PROGRAM linearised_BGK_fix_point
use omp_lib
implicit none
!***********************************************************************
!   Constant and variable declarations
!***********************************************************************

double precision, parameter :: accom=1.d0 ! accommodation coefficient
double precision, parameter :: porosity=0.75d0
double precision, parameter :: eps=1.d-8
!=======================================================================
!       Constant and indexed, temporary variables
!=======================================================================
double precision,parameter :: PI=datan(1.d0)*4.d0
!double precision,parameter :: PI=3.1416d0
integer :: i,j,k    !index variable for physical space
integer :: l,m,n  !index variable for molecular velocity space
logical :: fileexist
!=======================================================================
!       Physical space configurations
!=======================================================================
integer, parameter :: fluidlayer=4
!integer, parameter :: Ny = 101 !number of point half height channel, including two end
!integer, parameter :: Nx = (Ny-1)*2+1  !length of the channel
!integer, parameter :: Nz = Ny
integer, parameter :: Ny=30
integer, parameter :: Nx=Ny+2*fluidlayer
integer, parameter :: Nz=Ny
integer, parameter :: translate=0
integer, parameter :: ghostLayer=2 !number of ghost layer at each boundary
!integer, parameter :: obstR=(Nx-1)/4     !2D square cylinder with porosity=0.75
!integer, parameter :: obstR = floor((Nx-1)*dsqrt((1.d0-porosity)/PI))     !2D circular cylinder
!integer, parameter :: obstR = floor( dsqrt(1-porosity)/2.0*(Nx-1) )       !2D square cylinder
integer, parameter :: obstR = floor((Nx-1)*((1.d0-porosity)*3.d0/4.d0/PI)**(1.d0/3.d0))     !3D sphere obstacle
!integer, parameter :: obstR = floor( (1-porosity)**(1.d0/3.d0)/2.0*(Nx-1) )       !3D cube obstacle
integer, parameter :: obstX = ghostLayer+Ny-2   !Check
integer, parameter :: obstY = ghostLayer-1      !Check
integer, parameter :: obstZ = ghostLayer-1      !Check
integer, parameter :: Nxmin=1-ghostLayer,Nxmax=Nx+ghostLayer
integer, parameter :: Nymin=1-ghostLayer,Nymax=Ny+ghostLayer
integer, parameter :: Nzmin=1-ghostLayer,Nzmax=Nz+ghostLayer
integer, parameter :: Nxtotal=Nxmax-Nxmin+1
integer, parameter :: Nytotal=Nymax-Nymin+1
integer, parameter :: Nztotal=Nzmax-Nzmin+1
integer, parameter :: Ntotal=Nxtotal*Nytotal*Nztotal
integer, parameter :: Nxytotal=Nxtotal*Nytotal

integer, parameter:: fluid = 0, solid = 1, ghost=4
integer, parameter:: wallE = 20, wallW = 21, wallN = 22, wallS = 23, wallF = 24, wallB = 25 !(1-side wall )
integer, parameter:: wallEN = 30, wallWN = 31, wallES = 32, wallWS = 33 !(2-side wall)
integer, parameter:: wallNF = 40, wallNB = 41, wallSF = 42, wallSB = 43 !(2-side wall)
integer, parameter:: wallEF = 50, wallEB = 51, wallWF = 52, wallWB = 53 !(2-side wall)
integer, parameter:: wallENF = 60, wallWNF = 61, wallWSF = 62, wallESF = 63 !(3-side wall)
integer, parameter:: wallENB = 70, wallWNB = 71, wallWSB = 72, wallESB = 73 !(3-side wall)
integer :: NneighborSolid, NneighborFluid, Nstencil, icount, iteration1
integer :: Nstencil1, Nstencil2, Nstencil3, Nstencil4 !(group 1-4)
integer :: Nstencil5, Nstencil6, Nstencil7, Nstencil8 !(group 5-8)
double precision :: RhoWall,RhoWall2,RhoWall3
!integer :: itemp
!integer, dimension(fluidlayer+1:Nx-fluidlayer) :: itemp
Character(kind=1)::ctemp
!=======================================================================
!       Flow configuration configurations
!=======================================================================

integer, parameter :: column=fluidlayer/2 ! layer to extract flow rate
double precision,parameter :: PressDrop=1.0d-3 !Pressure drop
integer, parameter :: iteration_max=10000000,iteration_min=10,interval=100
double precision :: ds = 1.d0/(Nx-1) !uniform grid spacing in physical space
!double precision :: dt  !uniform grid spacing in  time space
double precision :: Kn,mass,mass2,permeability
double precision :: DiffFlux, error, fEq, mu
double precision :: start_time, end_time

integer, parameter :: nKn=1! number of Kn case 36
integer :: iKn
!double precision, parameter, DIMENSION(1:nKn)::seriesKn=(/10.d-1,9.d-1,8.d-1,7.d-1,6.d-1,5.d-1,4.d-1,3.d-1,2.d-1 &
!                                                        & ,10.d-2,9.d-2,8.d-2,7.d-2,6.d-2,5.d-2,4.d-2,3.d-2,2.d-2 &
!                                                        & ,10.d-3,9.d-3,8.d-3,7.d-3,6.d-3,5.d-3,4.d-3,3.d-3,2.d-3 &
!                                                        & ,10.d-4,9.d-4,8.d-4,7.d-4,6.d-4,5.d-4,4.d-4,3.d-4,2.d-4 /)
double precision, parameter, DIMENSION(1:nKn)::seriesKn=(/1.d-1/)

!=======================================================================
!       Molecular velocity space configurations
!=======================================================================
integer, parameter :: Nc_fundamental=2 !Number of fundamental molecular velocity
integer, parameter :: Nc=(2*Nc_fundamental)**3 !Number of molecular velocity
!integer, parameter :: Nc4=Nc/4 !Quarter-number of molecular velocity
integer, parameter :: Nc8=Nc/8 !Quarter-number of molecular velocity
double precision, DIMENSION (1:Nc_fundamental) :: xi,weight1D !abscissae and weighting Hermite quadrature
double precision, DIMENSION (1:Nc8) :: cx, cy, cz, w !molecular velocity and weighting
!integer, DIMENSION (1:Nc2) :: oppositeX, oppositeY! specular wall's normal vector in X, Y, Z direction
!=======================================================================
!       Microscopic and macroscopic parameters
!=======================================================================

integer :: nWall, nCorner
integer, DIMENSION(:), ALLOCATABLE :: vecWall
integer, DIMENSION(:), ALLOCATABLE :: dir1, dir2, dir3, dir4
integer, DIMENSION(:), ALLOCATABLE :: dir5, dir6, dir7, dir8
double precision, DIMENSION(:,:), ALLOCATABLE :: coef1, coef2, coef3, coef4
double precision, DIMENSION(:,:), ALLOCATABLE :: coef5, coef6, coef7, coef8

integer, DIMENSION(:), ALLOCATABLE :: image, which_corner
double precision, DIMENSION(Nc8):: f1wZ,f2wZ,f3wZ,f4wZ,f5wZ,f6wZ,f7wZ,f8wZ
double precision, DIMENSION(:,:), ALLOCATABLE :: f1,f2,f3,f4,f5,f6,f7,f8
double precision, DIMENSION(:,:), ALLOCATABLE :: f1w,f2w,f3w,f4w,f5w,f6w,f7w,f8w
double precision, DIMENSION(:), ALLOCATABLE :: Rho, Ux, Uy, Uz
ALLOCATE(image(Ntotal),which_corner(Ntotal))
ALLOCATE(f1(Ntotal,Nc8),f2(Ntotal,Nc8),f3(Ntotal,Nc8),f4(Ntotal,Nc8))   !Group1-4
ALLOCATE(f5(Ntotal,Nc8),f6(Ntotal,Nc8),f7(Ntotal,Nc8),f8(Ntotal,Nc8))   !Group5-8
ALLOCATE(Rho(Ntotal), Ux(Ntotal), Uy(Ntotal), Uz(Ntotal))
!------------------------------------------------------------------------
!           Fundamental molecular velocity
!------------------------------------------------------------------------
xi(1) = dsqrt(3.d0-dsqrt(6.d0))           !fundamental abscissae
xi(2) = dsqrt(3.d0+dsqrt(6.d0))
weight1D(1) = (3.d0+dsqrt(6.d0))/12.d0    !fundamental weighting
weight1D(2) = (3.d0-dsqrt(6.d0))/12.d0
!------------------------------------------------------------------------
!           Two-dimensional Hermite quadrature: tensor production formulae
!           Mapping from  2D array to 1D array
!------------------------------------------------------------------------
Do l=1,Nc_fundamental
    Do m=1,Nc_fundamental
        Do n=1,Nc_fundamental
            k=n+(m-1)*Nc_fundamental+(l-1)*Nc_fundamental**2
            cx(k) = xi(n)
            cy(k) = xi(m)
            cz(k) = xi(l)
            w(k) = weight1D(n)*weight1D(m)*weight1D(l)
        End do
    End do
End do
!***********************************************************************
!   Initialisation
!***********************************************************************
    inquire(file="Results.dat",exist=fileexist)
    if (.not.fileexist) then
        open(22,file="Results.dat",STATUS='NEW')
        write(22,*) ' TITLE=" Results"'
        write(22,*) ' VARIABLES=Kn,mass,permeability,error,iteration,wall_time_hours'
        write(22,*) ' ZONE T="final"'
        close(22)
    end if
!=======================================================================
!       Physical image of flow+solid+wall domain
!=======================================================================
image = fluid
nWall=0
which_corner=0
nCorner=0
!------------------------------------------------------------------------
!           Prefill ghosts layer by ghost
!------------------------------------------------------------------------
If (ghostLayer>0) then
Do k=Nzmin,Nzmax
    Do j=Nymin,Nymax
        Do i=Nxmin,Nxmax
            If ((i<1).OR.(i>Nx).OR.(j<1).OR.(j>Ny).OR.(k<1).OR.(k>Nz)) then
                l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
                image(l) = ghost
            End if
        Enddo
    Enddo
Enddo    
End if
!------------------------------------------------------------------------
!           Read from tomography file  !!!Check!!!
!------------------------------------------------------------------------
!Open(200,file='Gambier300^3.raw',status='OLD',form='unformatted',ACCESS="STREAM")
!Open(200,file='SandPack300^3.raw',status='OLD',form='unformatted',ACCESS="STREAM")
!Open(200,file='BeadPack300^3.raw',status='OLD',form='unformatted',ACCESS="STREAM")
!Open(200,file='Castlegate300^3.raw',status='OLD',form='unformatted',ACCESS="STREAM")
!Open(200,file='shale400^3.raw',status='OLD',form='unformatted',ACCESS="STREAM")


!Do k=1,Ny
!    Do j=1,Ny
!!Do k=fluidlayer+1,300+fluidlayer 
!!   Do j=fluidlayer+1,300+fluidlayer 
!        Do i=fluidlayer+1,Ny+fluidlayer  
!            read(200) ctemp       
!!           l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
!!       if ((k<=Nz).AND.(j<=Ny).AND.(i<=Nx-fluidlayer)) then
!            l=(i-translate+ghostLayer)+(j-translate+ghostLayer-1)*Nxtotal+(k-translate+ghostLayer-1)*Nxytotal
!        
!            if ((k>=1+translate).AND.(k<=Nz+translate).AND.(j>=1+translate).AND.(j<=Ny+translate).AND.&
!                 & (i>=fluidlayer+1+translate).AND.(i<=Nx-fluidlayer+translate)) then
!                if (ichar(ctemp)==0) then
!                    image(l) = fluid
!                else if (ichar(ctemp)==1) then
!                    image(l) = solid
!                else 
!                    write(*,*) "invalid image data file"
!                    stop 0
!                endif
!            endif
!        Enddo  
!    Enddo
!Enddo
!Close(200)


! DEBUG
!image = fluid
do k=11,20
    do j=11,20
        do i=15,28
            l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
            image(l) = solid
        Enddo
    Enddo
Enddo

!------------------------------------------------------------------------
!           Drill to obtain at least two layer of fluid in pore
!------------------------------------------------------------------------

Do k=1,Nz
  Do j=1,Ny
    Do i=2,Nx-1
    l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
    If (image(l)==fluid) then
      If (((image(l-1)/=fluid).AND.(image(l+1)/=fluid)).OR.((image(l-Nxtotal)/=fluid).AND.(image(l+Nxtotal)/=fluid))&
      & .OR.((image(l-Nxytotal)/=fluid).AND.(image(l+Nxytotal)/=fluid))) then
        If (j==Ny) then
          If (k==Nz) then
            image(l+1)=fluid
            image(l-Nxtotal)=fluid
            image(l+1-Nxtotal)=fluid
            image(l-Nxytotal)=fluid
            image(l+1-Nxytotal)=fluid
            image(l-Nxtotal-Nxytotal)=fluid
            image(l+1-Nxtotal-Nxytotal)=fluid  
          else
            image(l+1)=fluid
            image(l-Nxtotal)=fluid
            image(l+1-Nxtotal)=fluid
            image(l+Nxytotal)=fluid
            image(l+1+Nxytotal)=fluid
            image(l-Nxtotal+Nxytotal)=fluid
            image(l+1-Nxtotal+Nxytotal)=fluid     
          Endif 
        else if (k==Nz)  then
            image(l+1)=fluid
            image(l+Nxtotal)=fluid
            image(l+1+Nxtotal)=fluid
            image(l-Nxytotal)=fluid
            image(l+1-Nxytotal)=fluid
            image(l+Nxtotal-Nxytotal)=fluid
            image(l+1+Nxtotal-Nxytotal)=fluid       
        else
            image(l+1)=fluid
            image(l+Nxtotal)=fluid
            image(l+1+Nxtotal)=fluid
            image(l+Nxytotal)=fluid
            image(l+1+Nxytotal)=fluid
            image(l+Nxtotal+Nxytotal)=fluid
            image(l+1+Nxtotal+Nxytotal)=fluid
        endif           
      Endif     
    End if
    Enddo
  Enddo
Enddo   

!------------------------------------------------------------------------
!           Prefill ghosts layer by fluid
!------------------------------------------------------------------------
If (ghostLayer>0) then
Do k=Nzmin,Nzmax
    Do j=Nymin,Nymax
        Do i=Nxmin,Nxmax
            If ((i<1).OR.(i>Nx).OR.(j<1).OR.(j>Ny).OR.(k<1).OR.(k>Nz)) then
                l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
                image(l) = fluid
            End if
        Enddo
    Enddo
Enddo    
End if 

Do k=1,Nz
  Do j=1,Ny
    Do i=1,Nx
      l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
      If (image(l)==solid) then
    If (((image(l+1)==fluid).AND.(image(l-1)==fluid)).OR.((image(l+Nxtotal)==fluid).AND.(image(l-Nxtotal)==fluid))&
    & .OR.((image(l+Nxytotal)==fluid).AND.(image(l-Nxytotal)==fluid))) then ! 1-layer-thickness wall 
      image(l)=fluid
    Endif       
      Endif
    Enddo
  Enddo
Enddo  
! 

!------------------------------------------------------------------------
!           Prefill ghosts layer by ghost
!------------------------------------------------------------------------
If (ghostLayer>0) then
Do k=Nzmin,Nzmax
    Do j=Nymin,Nymax
        Do i=Nxmin,Nxmax
            If ((i<1).OR.(i>Nx).OR.(j<1).OR.(j>Ny).OR.(k<1).OR.(k>Nz)) then
                l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
                image(l) = ghost
            End if
        Enddo
    Enddo
Enddo    
End if

Do k=1,Nz
  Do j=1,Ny
    Do i=2,Nx-1
    l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
    If (image(l)==fluid) then
      If (((image(l-1)/=fluid).AND.(image(l+1)/=fluid)).OR.((image(l-Nxtotal)/=fluid).AND.(image(l+Nxtotal)/=fluid))&
      & .OR.((image(l-Nxytotal)/=fluid).AND.(image(l+Nxytotal)/=fluid))) then
    If (j==Ny) then
      If (k==Nz) then
        image(l+1)=fluid
        image(l-Nxtotal)=fluid
        image(l+1-Nxtotal)=fluid
        image(l-Nxytotal)=fluid
        image(l+1-Nxytotal)=fluid
        image(l-Nxtotal-Nxytotal)=fluid
        image(l+1-Nxtotal-Nxytotal)=fluid  
      else
        image(l+1)=fluid
        image(l-Nxtotal)=fluid
        image(l+1-Nxtotal)=fluid
        image(l+Nxytotal)=fluid
        image(l+1+Nxytotal)=fluid
        image(l-Nxtotal+Nxytotal)=fluid
        image(l+1-Nxtotal+Nxytotal)=fluid     
      Endif 
    else if (k==Nz)  then
        image(l+1)=fluid
        image(l+Nxtotal)=fluid
        image(l+1+Nxtotal)=fluid
        image(l-Nxytotal)=fluid
        image(l+1-Nxytotal)=fluid
        image(l+Nxtotal-Nxytotal)=fluid
        image(l+1+Nxtotal-Nxytotal)=fluid       
    else
        image(l+1)=fluid
        image(l+Nxtotal)=fluid
        image(l+1+Nxtotal)=fluid
        image(l+Nxytotal)=fluid
        image(l+1+Nxytotal)=fluid
        image(l+Nxtotal+Nxytotal)=fluid
        image(l+1+Nxtotal+Nxytotal)=fluid
    endif           
      Endif     
    End if
    Enddo
  Enddo
Enddo   

!------------------------------------------------------------------------
!           Prefill ghosts layer by fluid
!------------------------------------------------------------------------
If (ghostLayer>0) then
Do k=Nzmin,Nzmax
    Do j=Nymin,Nymax
        Do i=Nxmin,Nxmax
            If ((i<1).OR.(i>Nx).OR.(j<1).OR.(j>Ny).OR.(k<1).OR.(k>Nz)) then
                l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
                image(l) = fluid
            End if
        Enddo
    Enddo
Enddo    
End if

Do k=1,Nz
  Do j=1,Ny
    Do i=1,Nx
      l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
      If (image(l)==solid) then
    If (((image(l+1)==fluid).AND.(image(l-1)==fluid)).OR.((image(l+Nxtotal)==fluid).AND.(image(l-Nxtotal)==fluid))&
    & .OR.((image(l+Nxytotal)==fluid).AND.(image(l-Nxytotal)==fluid))) then ! 1-layer-thickness wall 
      image(l)=fluid
    Endif       
      Endif
    Enddo
  Enddo
Enddo  

!------------------------------------------------------------------------
!           Prefill ghosts layer by ghost
!------------------------------------------------------------------------
If (ghostLayer>0) then
Do k=Nzmin,Nzmax
    Do j=Nymin,Nymax
        Do i=Nxmin,Nxmax
            If ((i<1).OR.(i>Nx).OR.(j<1).OR.(j>Ny).OR.(k<1).OR.(k>Nz)) then
                l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
                image(l) = ghost
            End if
        Enddo
    Enddo
Enddo    
End if
!------------------------------------------------------------------------
!           Assign automatically the wall layer between fluid and solid 
!------------------------------------------------------------------------
Do k=1,Nz
Do j=1,Ny
    Do i=1,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==solid) then
            NneighborFluid=1 !(1-2-3-4 in D2Q9 corresponding to 2-3-5-7)
            If (image(l+1)==fluid)  NneighborFluid=NneighborFluid*2   !Check neighbor on the East(2)
            If (image(l-1)==fluid)  NneighborFluid=NneighborFluid*3   !Check neighbor on the West(3)
            If (image(l+Nxtotal)==fluid)  NneighborFluid=NneighborFluid*5   !Check neighbor on the North(5)
            If (image(l-Nxtotal)==fluid)  NneighborFluid=NneighborFluid*7   !Check neighbor on the South(7)
            If (image(l+Nxytotal)==fluid)  NneighborFluid=NneighborFluid*9   !Check neighbor on the Front(9)
            If (image(l-Nxytotal)==fluid)  NneighborFluid=NneighborFluid*11   !Check neighbor on the Back(11)

            SELECT Case (NneighborFluid)
                ! 1 fluid-neighbor wall
                CASE (2)
                    image(l)=wallE
                    nWall=nWall+1
                CASE (3)
                    image(l)=wallW
                    nWall=nWall+1
                CASE (5)
                    image(l)=wallN
                    nWall=nWall+1
                CASE (7)
                    image(l)=wallS
                    nWall=nWall+1
                CASE (9)
                    image(l)=wallF
                    nWall=nWall+1
                CASE (11)
                    image(l)=wallB
                    nWall=nWall+1
                ! 2 fluid-neighbor wall
                CASE (10)
                    image(l)=wallEN
                    nWall=nWall+1
                CASE (15)
                    image(l)=wallWN
                    nWall=nWall+1
                CASE (14)
                    image(l)=wallES
                    nWall=nWall+1
                CASE (21)
                    image(l)=wallWS
                    nWall=nWall+1
                CASE (18)
                    image(l)=wallEF
                    nWall=nWall+1
                CASE (22)
                    image(l)=wallEB
                    nWall=nWall+1
                CASE (27)
                    image(l)=wallWF
                    nWall=nWall+1
                CASE (33)
                    image(l)=wallWB
                    nWall=nWall+1
                CASE (45)
                    image(l)=wallNF
                    nWall=nWall+1
                CASE (55)
                    image(l)=wallNB
                    nWall=nWall+1
                CASE (63)
                    image(l)=wallSF
                    nWall=nWall+1
                CASE (77)
                    image(l)=wallSB
                    nWall=nWall+1
                ! 3 fluid-neighbor wall
                CASE (90)
                    image(l)=wallENF
                    nWall=nWall+1
                CASE (110)
                    image(l)=wallENB
                    nWall=nWall+1
                CASE (135)
                    image(l)=wallWNF
                    nWall=nWall+1
                CASE (165)
                    image(l)=wallWNB
                    nWall=nWall+1
                CASE (126)
                    image(l)=wallESF
                    nWall=nWall+1
                CASE (154)
                    image(l)=wallESB
                    nWall=nWall+1
                CASE (189)
                    image(l)=wallWSF
                    nWall=nWall+1
                CASE (231)
                    image(l)=wallWSB
                    nWall=nWall+1
            END SELECT
        Endif
    Enddo
Enddo
Enddo
!------------------------------------------------------------------------
!           Assign manually the wall layer  
!------------------------------------------------------------------------
!    Do i=1,Nx  
!        !bottom plane
!       j=1 
!       k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal
!        image(k)=wallN
!        nWall=nWall+1
!        !top plane 
!       j=Ny 
!       k=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal 
!        image(k)=wallS
!        nWall=nWall+1
!    End do
!------------------------------------------------------------------------
!           Create wall-type vectors
!------------------------------------------------------------------------
ALLOCATE(vecWall(nWall))
nWall=0
Do k=1,Nz
Do j=1,Ny
    Do i=1,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        SELECT Case (image(l))
            ! 1 fluid-neighbor wall
            CASE (wallE)
                nWall=nWall+1
                vecWall(nWall)=l
            CASE (wallW)
                nWall=nWall+1
                vecWall(nWall)=l
            CASE (wallN)
                nWall=nWall+1
                vecWall(nWall)=l
            CASE (wallS)
                nWall=nWall+1
                vecWall(nWall)=l
            CASE (wallF)
                nWall=nWall+1
                vecWall(nWall)=l
            CASE (wallB)
                nWall=nWall+1
                vecWall(nWall)=l
            ! 2 fluid-neighbor wall
            CASE (wallEN)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallWN)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallES)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallWS)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallEF)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallEB)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallWF)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallWB)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallNF)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallNB)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallSF)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallSB)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            ! 3 fluid-neighbor wall
            CASE (wallENF)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallENB)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallWNF)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallWNB)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallESF)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallESB)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallWSF)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
            CASE (wallWSB)
                nWall=nWall+1
                vecWall(nWall)=l
                nCorner=nCorner+1
                which_corner(l)=nCorner
        END SELECT
    Enddo
Enddo
Enddo
ALLOCATE(f1w(nCorner,Nc8),f2w(nCorner,Nc8),f3w(nCorner,Nc8),f4w(nCorner,Nc8) &
& ,f5w(nCorner,Nc8),f6w(nCorner,Nc8),f7w(nCorner,Nc8),f8w(nCorner,Nc8))
!------------------------------------------------------------------------
!           Assign the ghost layers
!------------------------------------------------------------------------
If (ghostLayer>0) then
Do k=Nzmin,Nzmax
    Do j=Nymin,Nymax
        Do i=Nxmin,Nxmax
            If ((i<1).OR.(i>Nx).OR.(j<1).OR.(j>Ny).OR.(k<1).OR.(k>Nz)) then
                l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
                image(l) = ghost
            End if
        Enddo
    Enddo
Enddo
End if

!=======================================================
!       Create advancing-point vectors
!=======================================================================

icount=0
Do k=2,Nz
Do j=2,Ny
    Do i=2,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
        End if
    End do
End do
End do
Nstencil1=icount !Number of fluid points in the numerical stencil

icount=0
Do k=2,Nz
Do j=2,Ny
    Do i=Nx-1,1,-1
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
        End if
    End do
End do
End do
Nstencil2=icount !Number of fluid points in the numerical stencil

icount=0
Do k=2,Nz
Do j=Ny-1,1,-1
    Do i=Nx-1,1,-1
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
        End if
    End do
End do
End do
Nstencil3=icount !Number of fluid points in the numerical stencil

icount=0
Do k=2,Nz
Do j=Ny-1,1,-1
    Do i=2,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
        End if
    End do
End do
End do
Nstencil4=icount !Number of fluid points in the numerical stencil

icount=0
Do k=Nz-1,1,-1
Do j=2,Ny
    Do i=2,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
        End if
    End do
End do
End do
Nstencil5=icount !Number of fluid points in the numerical stencil

icount=0
Do k=Nz-1,1,-1
Do j=2,Ny
    Do i=Nx-1,1,-1
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
        End if
    End do
End do
End do
Nstencil6=icount !Number of fluid points in the numerical stencil

icount=0
Do k=Nz-1,1,-1
Do j=Ny-1,1,-1
    Do i=Nx-1,1,-1
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
        End if
    End do
End do
End do
Nstencil7=icount !Number of fluid points in the numerical stencil

icount=0
Do k=Nz-1,1,-1
Do j=Ny-1,1,-1
    Do i=2,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
        End if
    End do
End do
End do
Nstencil8=icount !Number of fluid points in the numerical stencil

ALLOCATE(dir1(Nstencil1),dir2(Nstencil2),dir3(Nstencil3),dir4(Nstencil4))
ALLOCATE(dir5(Nstencil5),dir6(Nstencil6),dir7(Nstencil7),dir8(Nstencil8))

icount=0
Do k=2,Nz
Do j=2,Ny
    Do i=2,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
            dir1(icount)=l
        End if
    End do
End do
End do

icount=0
Do k=2,Nz
Do j=2,Ny
    Do i=Nx-1,1,-1
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
            dir2(icount)=l
        End if
    End do
End do
End do


icount=0
Do k=2,Nz
Do j=Ny-1,1,-1
    Do i=Nx-1,1,-1
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
            dir3(icount)=l
        End if
    End do
End do
End do


icount=0
Do k=2,Nz
Do j=Ny-1,1,-1
    Do i=2,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
            dir4(icount)=l
        End if
    End do
End do
End do


icount=0
Do k=Nz-1,1,-1
Do j=2,Ny
    Do i=2,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
            dir5(icount)=l
        End if
    End do
End do
End do


icount=0
Do k=Nz-1,1,-1
Do j=2,Ny
    Do i=Nx-1,1,-1
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
            dir6(icount)=l
        End if
    End do
End do
End do


icount=0
Do k=Nz-1,1,-1
Do j=Ny-1,1,-1
    Do i=Nx-1,1,-1
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
            dir7(icount)=l
        End if
    End do
End do
End do


icount=0
Do k=Nz-1,1,-1
Do j=Ny-1,1,-1
    Do i=2,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
            icount=icount+1
            dir8(icount)=l
        End if
    End do
End do
End do
!=======================================================================
!       Create matrix of coefficients used in numerical stencil
!       coefI(point, weight_x0n & weight_x1n & weight_x2n, weight_y0n & weight_y1n & weight_y2n)
!       coefII(point, weight_x0p & weight_x1p & weight_x2p, weight_y0n & weight_y1n & weight_y2n)
!       coefIII(point, weight_x0p & weight_x1p & weight_x2p, weight_y0p & weight_y1p & weight_y2p)
!       coefIV(point, weight_x0n & weight_x1n & weight_x2n, weight_y0p & weight_y1p & weight_y2p)
!       coefI(i,1)-coefI(i,2)-coefI(i,3)=0, coefI(i,4)-coefI(i,5)-coefI(i,6)=0
!=======================================================================
ALLOCATE(coef1(Nstencil1,9),coef2(Nstencil2,9),coef3(Nstencil3,9),coef4(Nstencil4,9))
ALLOCATE(coef5(Nstencil5,9),coef6(Nstencil6,9),coef7(Nstencil7,9),coef8(Nstencil8,9))
Do i=1,Nstencil1
    k=dir1(i)
    !2nd order of accuracy
    coef1(i,1)=1.5d0/ds !x0n
    coef1(i,2)=2.d0/ds  !x1n
    coef1(i,3)=-0.5d0/ds   !x2n
    coef1(i,4)=1.5d0/ds !y0n
    coef1(i,5)=2.d0/ds  !y1n
    coef1(i,6)=-0.5d0/ds   !y2n
    coef1(i,7)=1.5d0/ds !z0n
    coef1(i,8)=2.d0/ds  !z1n
    coef1(i,9)=-0.5d0/ds   !z2n
    !1st order of accuracy in x
    if ((image(k-2)==ghost).OR.(image(k-1)==wallE) &
    &.OR.(image(k-1)==wallES).OR.(image(k-1)==wallEN).OR.(image(k-1)==wallEF).OR.(image(k-1)==wallEB) &
    &.OR.(image(k-1)==wallENF).OR.(image(k-1)==wallENB).OR.(image(k-1)==wallESF).OR.(image(k-1)==wallESB)) then
        coef1(i,1)=1.d0/ds  !x0n
        coef1(i,2)=1.d0/ds  !x1n
        coef1(i,3)=0.d0     !x2n
    end if
    !1st order of accuracy in y
    if ((image(k-2*Nxtotal)==ghost).OR.(image(k-Nxtotal)==wallN) &
    &.OR.(image(k-Nxtotal)==wallEN).OR.(image(k-Nxtotal)==wallWN).OR.(image(k-Nxtotal)==wallNF).OR.(image(k-Nxtotal)==wallNB) &
    &.OR.(image(k-Nxtotal)==wallENF).OR.(image(k-Nxtotal)==wallENB).OR.(image(k-Nxtotal)==wallWNF).OR.(image(k-Nxtotal)==wallWNB)) then
        coef1(i,4)=1.d0/ds  !y0n
        coef1(i,5)=1.d0/ds  !y1n
        coef1(i,6)=0.d0     !y2n
    end if
    !1st order of accuracy in z
    if ((image(k-2*Nxytotal)==ghost).OR.(image(k-Nxytotal)==wallF) &
    &.OR.(image(k-Nxytotal)==wallNF).OR.(image(k-Nxytotal)==wallSF).OR.(image(k-Nxytotal)==wallEF).OR.(image(k-Nxytotal)==wallWF) &
    &.OR.(image(k-Nxytotal)==wallENF).OR.(image(k-Nxytotal)==wallWNF).OR.(image(k-Nxytotal)==wallESF) &
    &.OR.(image(k-Nxytotal)==wallWSF)) then
        coef1(i,7)=1.d0/ds  !z0n
        coef1(i,8)=1.d0/ds  !z1n
        coef1(i,9)=0.d0     !z2n
    end if
End do

Do i=1,Nstencil2
    k=dir2(i)
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
    if ((image(k+2)==ghost).OR.(image(k+1)==wallW) &
    &.OR.(image(k+1)==wallWS).OR.(image(k+1)==wallWN).OR.(image(k+1)==wallWF).OR.(image(k+1)==wallWB) &
    &.OR.(image(k+1)==wallWNF).OR.(image(k+1)==wallWNB).OR.(image(k+1)==wallWSF).OR.(image(k+1)==wallWSB)) then
        coef2(i,1)=-1.d0/ds  !x0n
        coef2(i,2)=-1.d0/ds  !x1n
        coef2(i,3)=0.d0     !x2n
    end if
    !1st order of accuracy in y
    if ((image(k-2*Nxtotal)==ghost).OR.(image(k-Nxtotal)==wallN) &
    &.OR.(image(k-Nxtotal)==wallEN).OR.(image(k-Nxtotal)==wallWN).OR.(image(k-Nxtotal)==wallNF).OR.(image(k-Nxtotal)==wallNB) &
    &.OR.(image(k-Nxtotal)==wallENF).OR.(image(k-Nxtotal)==wallENB).OR.(image(k-Nxtotal)==wallWNF) &
    &.OR.(image(k-Nxtotal)==wallWNB)) then
        coef2(i,4)=1.d0/ds  !y0n
        coef2(i,5)=1.d0/ds  !y1n
        coef2(i,6)=0.d0     !y2n
    end if
    !1st order of accuracy in z
    if ((image(k-2*Nxytotal)==ghost).OR.(image(k-Nxytotal)==wallF) &
    &.OR.(image(k-Nxytotal)==wallNF).OR.(image(k-Nxytotal)==wallSF).OR.(image(k-Nxytotal)==wallEF).OR.(image(k-Nxytotal)==wallWF) &
    &.OR.(image(k-Nxytotal)==wallENF).OR.(image(k-Nxytotal)==wallWNF).OR.(image(k-Nxytotal)==wallESF) &
    &.OR.(image(k-Nxytotal)==wallWSF)) then
        coef2(i,7)=1.d0/ds  !z0n
        coef2(i,8)=1.d0/ds  !z1n
        coef2(i,9)=0.d0     !z2n
    end if
End do

Do i=1,Nstencil3
    k=dir3(i)
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
    if ((image(k+2)==ghost).OR.(image(k+1)==wallW) &
    &.OR.(image(k+1)==wallWS).OR.(image(k+1)==wallWN).OR.(image(k+1)==wallWF).OR.(image(k+1)==wallWB) &
    &.OR.(image(k+1)==wallWNF).OR.(image(k+1)==wallWNB).OR.(image(k+1)==wallWSF).OR.(image(k+1)==wallWSB)) then
        coef3(i,1)=-1.d0/ds  !x0n
        coef3(i,2)=-1.d0/ds  !x1n
        coef3(i,3)=0.d0     !x2n
    end if
    !1st order of accuracy in y
    if ((image(k+2*Nxtotal)==ghost).OR.(image(k+Nxtotal)==wallS) &
    &.OR.(image(k+Nxtotal)==wallES).OR.(image(k+Nxtotal)==wallWS).OR.(image(k+Nxtotal)==wallSF).OR.(image(k+Nxtotal)==wallSB) &
    &.OR.(image(k+Nxtotal)==wallESF).OR.(image(k+Nxtotal)==wallESB).OR.(image(k+Nxtotal)==wallWSF).OR.(image(k+Nxtotal)==wallWSB)) then
        coef3(i,4)=-1.d0/ds  !y0n
        coef3(i,5)=-1.d0/ds  !y1n
        coef3(i,6)=0.d0     !y2n
    end if
    !1st order of accuracy in z
    if ((image(k-2*Nxytotal)==ghost).OR.(image(k-Nxytotal)==wallF) &
    &.OR.(image(k-Nxytotal)==wallNF).OR.(image(k-Nxytotal)==wallSF).OR.(image(k-Nxytotal)==wallEF).OR.(image(k-Nxytotal)==wallWF) &
    &.OR.(image(k-Nxytotal)==wallENF).OR.(image(k-Nxytotal)==wallWNF).OR.(image(k-Nxytotal)==wallESF) &
    & .OR.(image(k-Nxytotal)==wallWSF)) then
        coef3(i,7)=1.d0/ds  !z0n
        coef3(i,8)=1.d0/ds  !z1n
        coef3(i,9)=0.d0     !z2n
    end if
End do

Do i=1,Nstencil4
    k=dir4(i)
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
    if ((image(k-2)==ghost).OR.(image(k-1)==wallE) &
    &.OR.(image(k-1)==wallES).OR.(image(k-1)==wallEN).OR.(image(k-1)==wallEF).OR.(image(k-1)==wallEB) &
    &.OR.(image(k-1)==wallENF).OR.(image(k-1)==wallENB).OR.(image(k-1)==wallESF).OR.(image(k-1)==wallESB)) then
        coef4(i,1)=1.d0/ds  !x0n
        coef4(i,2)=1.d0/ds  !x1n
        coef4(i,3)=0.d0     !x2n
    end if
    !1st order of accuracy in y
    if ((image(k+2*Nxtotal)==ghost).OR.(image(k+Nxtotal)==wallS) &
    &.OR.(image(k+Nxtotal)==wallES).OR.(image(k+Nxtotal)==wallWS).OR.(image(k+Nxtotal)==wallSF).OR.(image(k+Nxtotal)==wallSB) &
    &.OR.(image(k+Nxtotal)==wallESF).OR.(image(k+Nxtotal)==wallESB).OR.(image(k+Nxtotal)==wallWSF) &
    & .OR.(image(k+Nxtotal)==wallWSB)) then
        coef4(i,4)=-1.d0/ds  !y0n
        coef4(i,5)=-1.d0/ds  !y1n
        coef4(i,6)=0.d0     !y2n
    end if
    !1st order of accuracy in z
    if ((image(k-2*Nxytotal)==ghost).OR.(image(k-Nxytotal)==wallF) &
    &.OR.(image(k-Nxytotal)==wallNF).OR.(image(k-Nxytotal)==wallSF).OR.(image(k-Nxytotal)==wallEF).OR.(image(k-Nxytotal)==wallWF) &
    &.OR.(image(k-Nxytotal)==wallENF).OR.(image(k-Nxytotal)==wallWNF).OR.(image(k-Nxytotal)==wallESF) &
    & .OR.(image(k-Nxytotal)==wallWSF)) then
        coef4(i,7)=1.d0/ds  !z0n
        coef4(i,8)=1.d0/ds  !z1n
        coef4(i,9)=0.d0     !z2n
    end if
End do

Do i=1,Nstencil5
    k=dir5(i)
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
    if ((image(k-2)==ghost).OR.(image(k-1)==wallE) &
    &.OR.(image(k-1)==wallES).OR.(image(k-1)==wallEN).OR.(image(k-1)==wallEF).OR.(image(k-1)==wallEB) &
    &.OR.(image(k-1)==wallENF).OR.(image(k-1)==wallENB).OR.(image(k-1)==wallESF).OR.(image(k-1)==wallESB)) then
        coef5(i,1)=1.d0/ds  !x0n
        coef5(i,2)=1.d0/ds  !x1n
        coef5(i,3)=0.d0     !x2n
    end if
    !1st order of accuracy in y
    if ((image(k-2*Nxtotal)==ghost).OR.(image(k-Nxtotal)==wallN) &
    &.OR.(image(k-Nxtotal)==wallEN).OR.(image(k-Nxtotal)==wallWN).OR.(image(k-Nxtotal)==wallNF).OR.(image(k-Nxtotal)==wallNB) &
    &.OR.(image(k-Nxtotal)==wallENF).OR.(image(k-Nxtotal)==wallENB).OR.(image(k-Nxtotal)==wallWNF) &
    & .OR.(image(k-Nxtotal)==wallWNB)) then
        coef5(i,4)=1.d0/ds  !y0n
        coef5(i,5)=1.d0/ds  !y1n
        coef5(i,6)=0.d0     !y2n
    end if
    !1st order of accuracy in z
    if ((image(k+2*Nxytotal)==ghost).OR.(image(k+Nxytotal)==wallB) &
    &.OR.(image(k+Nxytotal)==wallNB).OR.(image(k+Nxytotal)==wallSB).OR.(image(k+Nxytotal)==wallEB).OR.(image(k+Nxytotal)==wallWB) &
    &.OR.(image(k+Nxytotal)==wallENB).OR.(image(k+Nxytotal)==wallWNB).OR.(image(k+Nxytotal)==wallESB) &
    & .OR.(image(k+Nxytotal)==wallWSB)) then
        coef5(i,7)=-1.d0/ds  !z0n
        coef5(i,8)=-1.d0/ds  !z1n
        coef5(i,9)=0.d0     !z2n
    end if
End do

Do i=1,Nstencil6
    k=dir6(i)
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
    if ((image(k+2)==ghost).OR.(image(k+1)==wallW) &
    &.OR.(image(k+1)==wallWS).OR.(image(k+1)==wallWN).OR.(image(k+1)==wallWF).OR.(image(k+1)==wallWB) &
    &.OR.(image(k+1)==wallWNF).OR.(image(k+1)==wallWNB).OR.(image(k+1)==wallWSF).OR.(image(k+1)==wallWSB)) then
        coef6(i,1)=-1.d0/ds  !x0n
        coef6(i,2)=-1.d0/ds  !x1n
        coef6(i,3)=0.d0     !x2n
    end if
    !1st order of accuracy in y
    if ((image(k-2*Nxtotal)==ghost).OR.(image(k-Nxtotal)==wallN) &
    &.OR.(image(k-Nxtotal)==wallEN).OR.(image(k-Nxtotal)==wallWN).OR.(image(k-Nxtotal)==wallNF).OR.(image(k-Nxtotal)==wallNB) &
    &.OR.(image(k-Nxtotal)==wallENF).OR.(image(k-Nxtotal)==wallENB).OR.(image(k-Nxtotal)==wallWNF) &
    & .OR.(image(k-Nxtotal)==wallWNB)) then
        coef6(i,4)=1.d0/ds  !y0n
        coef6(i,5)=1.d0/ds  !y1n
        coef6(i,6)=0.d0     !y2n
    end if
    !1st order of accuracy in z
    if ((image(k+2*Nxytotal)==ghost).OR.(image(k+Nxytotal)==wallB) &
    &.OR.(image(k+Nxytotal)==wallNB).OR.(image(k+Nxytotal)==wallSB).OR.(image(k+Nxytotal)==wallEB).OR.(image(k+Nxytotal)==wallWB) &
    &.OR.(image(k+Nxytotal)==wallENB).OR.(image(k+Nxytotal)==wallWNB).OR.(image(k+Nxytotal)==wallESB) &
    & .OR.(image(k+Nxytotal)==wallWSB)) then
        coef6(i,7)=-1.d0/ds  !z0n
        coef6(i,8)=-1.d0/ds  !z1n
        coef6(i,9)=0.d0     !z2n
    end if
End do

Do i=1,Nstencil7
    k=dir7(i)
    !2nd order of accuracy
    coef7(i,1)=-1.5d0/ds !x0n
    coef7(i,2)=-2.d0/ds  !x1n
    coef7(i,3)=0.5d0/ds   !x2n
    coef7(i,4)=-1.5d0/ds !y0n
    coef7(i,5)=-2.d0/ds  !y1n
    coef7(i,6)=0.5d0/ds   !y2n
    coef7(i,7)=-1.5d0/ds !z0n
    coef7(i,8)=-2.d0/ds  !z1n
    coef7(i,9)=0.5d0/ds   !z2n
    !1st order of accuracy in x
    if ((image(k+2)==ghost).OR.(image(k+1)==wallW) &
    &.OR.(image(k+1)==wallWS).OR.(image(k+1)==wallWN).OR.(image(k+1)==wallWF).OR.(image(k+1)==wallWB) &
    &.OR.(image(k+1)==wallWNF).OR.(image(k+1)==wallWNB).OR.(image(k+1)==wallWSF).OR.(image(k+1)==wallWSB)) then
        coef7(i,1)=-1.d0/ds  !x0n
        coef7(i,2)=-1.d0/ds  !x1n
        coef7(i,3)=0.d0     !x2n
    end if
    !1st order of accuracy in y
    if ((image(k+2*Nxtotal)==ghost).OR.(image(k+Nxtotal)==wallS) &
    &.OR.(image(k+Nxtotal)==wallES).OR.(image(k+Nxtotal)==wallWS).OR.(image(k+Nxtotal)==wallSF).OR.(image(k+Nxtotal)==wallSB) &
    &.OR.(image(k+Nxtotal)==wallESF).OR.(image(k+Nxtotal)==wallESB).OR.(image(k+Nxtotal)==wallWSF) &
    & .OR.(image(k+Nxtotal)==wallWSB)) then
        coef7(i,4)=-1.d0/ds  !y0n
        coef7(i,5)=-1.d0/ds  !y1n
        coef7(i,6)=0.d0     !y2n
    end if
    !1st order of accuracy in z
    if ((image(k+2*Nxytotal)==ghost).OR.(image(k+Nxytotal)==wallB) &
    &.OR.(image(k+Nxytotal)==wallNB).OR.(image(k+Nxytotal)==wallSB).OR.(image(k+Nxytotal)==wallEB).OR.(image(k+Nxytotal)==wallWB) &
    &.OR.(image(k+Nxytotal)==wallENB).OR.(image(k+Nxytotal)==wallWNB).OR.(image(k+Nxytotal)==wallESB) &
    & .OR.(image(k+Nxytotal)==wallWSB)) then
        coef7(i,7)=-1.d0/ds  !z0n
        coef7(i,8)=-1.d0/ds  !z1n
        coef7(i,9)=0.d0     !z2n
    end if
End do

Do i=1,Nstencil8
    k=dir8(i)
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
    if ((image(k-2)==ghost).OR.(image(k-1)==wallE) &
    &.OR.(image(k-1)==wallES).OR.(image(k-1)==wallEN).OR.(image(k-1)==wallEF).OR.(image(k-1)==wallEB) &
    &.OR.(image(k-1)==wallENF).OR.(image(k-1)==wallENB).OR.(image(k-1)==wallESF).OR.(image(k-1)==wallESB)) then
        coef8(i,1)=1.d0/ds  !x0n
        coef8(i,2)=1.d0/ds  !x1n
        coef8(i,3)=0.d0     !x2n
    end if
    !1st order of accuracy in y
    if ((image(k+2*Nxtotal)==ghost).OR.(image(k+Nxtotal)==wallS) &
    &.OR.(image(k+Nxtotal)==wallES).OR.(image(k+Nxtotal)==wallWS).OR.(image(k+Nxtotal)==wallSF).OR.(image(k+Nxtotal)==wallSB) &
    &.OR.(image(k+Nxtotal)==wallESF).OR.(image(k+Nxtotal)==wallESB).OR.(image(k+Nxtotal)==wallWSF) &
    & .OR.(image(k+Nxtotal)==wallWSB)) then
        coef8(i,4)=-1.d0/ds  !y0n
        coef8(i,5)=-1.d0/ds  !y1n
        coef8(i,6)=0.d0     !y2n
    end if
    !1st order of accuracy in z
    if ((image(k+2*Nxytotal)==ghost).OR.(image(k+Nxytotal)==wallB) &
    &.OR.(image(k+Nxytotal)==wallNB).OR.(image(k+Nxytotal)==wallSB).OR.(image(k+Nxytotal)==wallEB).OR.(image(k+Nxytotal)==wallWB) &
    &.OR.(image(k+Nxytotal)==wallENB).OR.(image(k+Nxytotal)==wallWNB).OR.(image(k+Nxytotal)==wallESB) &
    & .OR.(image(k+Nxytotal)==wallWSB)) then
        coef8(i,7)=-1.d0/ds  !z0n
        coef8(i,8)=-1.d0/ds  !z1n
        coef8(i,9)=0.d0     !z2n
    end if
End do
!=======================================================================
!     Constant diffusive reflected flux at wall boundary
!=======================================================================
    DiffFlux=0.d0
    Do l=1,Nc8
        DiffFlux=DiffFlux+cz(l)*w(l)
    Enddo
    DiffFlux=DiffFlux*4 !Check
!=======================================================================
!       Macroscopic and microscopic parameters
!=======================================================================
!Do iKn=1,nKn
Do iKn=2,2

Kn=seriesKn(iKn)*dsqrt(2.d0/PI)
!mu=sqrt(PI)/2/Kn
mu=1.d0/Kn
error=1.d0
mass=1.d0
iteration1=0


Rho=0.d0 ! For linearised fomula
Ux=0.d0
Uy=0.d0
Uz=0.d0
f1w=0.d0
f2w=0.d0
f3w=0.d0
f4w=0.d0
f5w=0.d0
f6w=0.d0
f7w=0.d0
f8w=0.d0
f1wZ=0.d0
f2wZ=0.d0
f3wZ=0.d0
f4wZ=0.d0
f5wZ=0.d0
f6wZ=0.d0
f7wZ=0.d0
f8wZ=0.d0
f1=0.d0
Do k=1,Nz
Do j=1,Ny
     Do i=1,Nx
         l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
         if (image(l)/=solid) then
!             Rho(l)=PressDrop*(i/2.d0-Nx)/Nx
        Rho(l)= 0.d0
         end if
         f1(l,:)=w(:)*Rho(l) !Check
     Enddo
Enddo
Enddo
f2=f1
f3=f1
f4=f1
f5=f1
f6=f1
f7=f1
f8=f1


!***********************************************************************
!   Print out the input parameter
!***********************************************************************
    open(10,file='Test.dat',STATUS="REPLACE")
        write(10,*) 'Velocity space D3Q64: l, cx, cy, cz, w'
        Do l=1,Nc8
            write(10,*) l, cx(l), cy(l), cz(l), w(l)
        Enddo
    close(10)


!***********************************************************************
!   Main procedure
!***********************************************************************
call cpu_time(start_time)

Do WHILE ((error>eps).AND.(iteration1<iteration_max))
!Do WHILE ((iteration1<iteration_min))
!Do WHILE ((error>eps))
!=======================================================================
!     Reset summational variables
!=======================================================================

!=======================================================================
!     Collision and streaming
!=======================================================================

!$OMP PARALLEL & 
!$OMP DEFAULT(NONE) &
!$OMP PRIVATE(i,j,k,l) &
!$OMP PRIVATE(fEq,RhoWall,RhoWall2,RhoWall3) &
!$OMP PRIVATE(f1wZ,f2wZ,f3wZ,f4wZ,f5wZ,f6wZ,f7wZ,f8wZ) &

!!$OMP PRIVATE() &
!$OMP SHARED(nWall,vecWall) &
!$OMP SHARED(cx,cy,cz,w,mu) &
!$OMP SHARED(Nstencil1,Nstencil2,Nstencil3,Nstencil4) &
!$OMP SHARED(Nstencil5,Nstencil6,Nstencil7,Nstencil8) &
!$OMP SHARED(dir1,dir2,dir3,dir4) &
!$OMP SHARED(dir5,dir6,dir7,dir8) &
!$OMP SHARED(coef1,coef2,coef3,coef4) &
!$OMP SHARED(coef5,coef6,coef7,coef8) &
!$OMP SHARED(image,which_corner) &
!$OMP SHARED(f1,f2,f3,f4,f5,f6,f7,f8) &
!$OMP SHARED(f1w,f2w,f3w,f4w,f5w,f6w,f7w,f8w) &

!$OMP SHARED(Rho,Ux,Uy,Uz,DiffFlux)

!!$OMP SHARED() &
!!$OMP SHARED()
!------------------------------------------------------------------------
!           In the 1st group of direction cx>0 & cy>0 & cz>0
!------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC) 
    Do l=1,Nc8
        Do i=1,Nstencil1
            k=dir1(i)
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction x-y
            If ((image(k-Nxtotal)==wallEN).OR.(image(k-Nxtotal)==wallENF).OR.(image(k-Nxtotal)==wallENB)) then
                f1(k-Nxtotal,l)=f1w(which_corner(k-Nxtotal),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(k-Nxytotal)==wallNF).OR.(image(k-Nxytotal)==wallEF).OR.(image(k-Nxytotal)==wallWNF) &
            & .OR.(image(k-Nxytotal)==wallESF)) then
                f1(k-Nxytotal,l)=f1w(which_corner(k-Nxytotal),l)
            else if (image(k-Nxytotal)==wallENF) then
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
            k=dir2(i)
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected x&y-direction
            If ((image(k-Nxtotal)==wallWN).OR.(image(k-Nxtotal)==wallWNF).OR.(image(k-Nxtotal)==wallWNB)) then
                f2(k-Nxtotal,l)=f2w(which_corner(k-Nxtotal),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(k-Nxytotal)==wallNF).OR.(image(k-Nxytotal)==wallWF).OR.(image(k-Nxytotal)==wallENF) &
            & .OR.(image(k-Nxytotal)==wallWSF)) then
                f2(k-Nxytotal,l)=f2w(which_corner(k-Nxytotal),l)
            else if (image(k-Nxytotal)==wallWNF) then
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
            k=dir3(i)
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction
            If ((image(k+Nxtotal)==wallWS).OR.(image(k+Nxtotal)==wallWSF).OR.(image(k+Nxtotal)==wallWSB)) then
                f3(k+Nxtotal,l)=f3w(which_corner(k+Nxtotal),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(k-Nxytotal)==wallSF).OR.(image(k-Nxytotal)==wallWF).OR.(image(k-Nxytotal)==wallWNF) &
            & .OR.(image(k-Nxytotal)==wallESF)) then
                f3(k-Nxytotal,l)=f3w(which_corner(k-Nxytotal),l)
            else if (image(k-Nxytotal)==wallWSF) then
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
            k=dir4(i)
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction
            If ((image(k+Nxtotal)==wallES).OR.(image(k+Nxtotal)==wallESF).OR.(image(k+Nxtotal)==wallESB)) then
                f4(k+Nxtotal,l)=f4w(which_corner(k+Nxtotal),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(k-Nxytotal)==wallSF).OR.(image(k-Nxytotal)==wallEF).OR.(image(k-Nxytotal)==wallENF) &
            & .OR.(image(k-Nxytotal)==wallWSF)) then
                f4(k-Nxytotal,l)=f4w(which_corner(k-Nxytotal),l)
            else if (image(k-Nxytotal)==wallESF) then
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
            k=dir5(i)
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction x-y
            If ((image(k-Nxtotal)==wallEN).OR.(image(k-Nxtotal)==wallENF).OR.(image(k-Nxtotal)==wallENB)) then
                f5(k-Nxtotal,l)=f5w(which_corner(k-Nxtotal),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(k+Nxytotal)==wallNB).OR.(image(k+Nxytotal)==wallEB).OR.(image(k+Nxytotal)==wallESB) &
            & .OR.(image(k+Nxytotal)==wallWNB)) then
                f5(k+Nxytotal,l)=f5w(which_corner(k+Nxytotal),l)
            else if (image(k+Nxytotal)==wallENB) then
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
            k=dir6(i)
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected x&y-direction
            If ((image(k-Nxtotal)==wallWN).OR.(image(k-Nxtotal)==wallWNF).OR.(image(k-Nxtotal)==wallWNB)) then
                f6(k-Nxtotal,l)=f6w(which_corner(k-Nxtotal),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(k+Nxytotal)==wallNB).OR.(image(k+Nxytotal)==wallWB).OR.(image(k+Nxytotal)==wallENB) &
            & .OR.(image(k+Nxytotal)==wallWSB)) then
                f6(k+Nxytotal,l)=f6w(which_corner(k+Nxytotal),l)
            else if (image(k+Nxytotal)==wallWNB) then
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
            k=dir7(i)
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction
            If ((image(k+Nxtotal)==wallWS).OR.(image(k+Nxtotal)==wallWSF).OR.(image(k+Nxtotal)==wallWSB)) then
                f7(k+Nxtotal,l)=f7w(which_corner(k+Nxtotal),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(k+Nxytotal)==wallSB).OR.(image(k+Nxytotal)==wallWB).OR.(image(k+Nxytotal)==wallWNB) &
            & .OR.(image(k+Nxytotal)==wallESB)) then
                f7(k+Nxytotal,l)=f7w(which_corner(k+Nxytotal),l)
            else if (image(k+Nxytotal)==wallWSB) then
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
            k=dir8(i)
            ! Switch from reflected in x-direction (default) to in y-direction
            ! only for group of velocity overlapped by two reflected direction
            If ((image(k+Nxtotal)==wallES).OR.(image(k+Nxtotal)==wallESF).OR.(image(k+Nxtotal)==wallESB)) then
                f8(k+Nxtotal,l)=f8w(which_corner(k+Nxtotal),l)
            End if
            ! Switch from reflected in x/y-direction to in z-direction
            ! only for group of velocity overlapped by two/three reflected x/y&z-direction
            If ((image(k+Nxytotal)==wallSB).OR.(image(k+Nxytotal)==wallEB).OR.(image(k+Nxytotal)==wallENB) &
            & .OR.(image(k+Nxytotal)==wallWSB)) then
                f8(k+Nxytotal,l)=f8w(which_corner(k+Nxytotal),l)
            else if (image(k+Nxytotal)==wallESB) then
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

!=======================================================================
!     Boundary condition on the flat wall
!=======================================================================
!$OMP DO SCHEDULE(STATIC)
    Do i=1,nWall
        k=vecWall(i)
        j=which_corner(k)
        RhoWall=0.d0
        RhoWall2=0.d0
        RhoWall3=0.d0
        SELECT CASE (image(k))
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
!=======================================================================
!     Boundary condition on inlet & outlet
!=======================================================================
!$OMP DO SCHEDULE(STATIC) 
Do k=1,Nz
    i=1
    Do j=1,Ny
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        !inlet
        f1(l,:)=f1(l+Nx-1,:)+w(:)*PressDrop
        f4(l,:)=f4(l+Nx-1,:)+w(:)*PressDrop
        f5(l,:)=f5(l+Nx-1,:)+w(:)*PressDrop
        f8(l,:)=f8(l+Nx-1,:)+w(:)*PressDrop
        !outlet
        f2(l+Nx-1,:)=f2(l,:)-w(:)*PressDrop
        f3(l+Nx-1,:)=f3(l,:)-w(:)*PressDrop
        f6(l+Nx-1,:)=f6(l,:)-w(:)*PressDrop
        f7(l+Nx-1,:)=f7(l,:)-w(:)*PressDrop
    End do
End do
!$OMP END DO    

! !$OMP DO SCHEDULE(STATIC) 
! Do k=1,Nz
!     i=1
!     Do j=1,Ny
!         l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
!         !inlet
!         f1(l,:)=w(:)*PressDrop
!         f4(l,:)=w(:)*PressDrop
!         f5(l,:)=w(:)*PressDrop
!         f8(l,:)=w(:)*PressDrop
!         !outlet
!         f2(l+Nx-1,:)=0.d0
!         f3(l+Nx-1,:)=0.d0
!         f6(l+Nx-1,:)=0.d0
!         f7(l+Nx-1,:)=0.d0
!     End do
! End do
! !$OMP END DO 
!=======================================================================
!     Boundary condition on symmetric planes y=
!=======================================================================
! !$OMP DO SCHEDULE(STATIC) 
! Do k=1,Nz
!     j=1 ! check if we applied on ghost layer or on the first fluid layer
!     Do i=1,Nx
!         l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
!         !bottom plane
!         f3(l,:)=2.d0*f3(l+Nxtotal,:)-f3(l+2*Nxtotal,:) !Check range is whole plane
!         f4(l,:)=2.d0*f4(l+Nxtotal,:)-f4(l+2*Nxtotal,:) !or only the gap-line for this
!         f7(l,:)=2.d0*f7(l+Nxtotal,:)-f7(l+2*Nxtotal,:) !extrapolation
!         f8(l,:)=2.d0*f8(l+Nxtotal,:)-f8(l+2*Nxtotal,:)
!         f2(l,:)=f3(l,:)
!         f1(l,:)=f4(l,:)
!         f6(l,:)=f7(l,:)
!         f5(l,:)=f8(l,:)
!     End do
! End do
! !$OMP END DO NOWAIT
! 
! !$OMP DO SCHEDULE(STATIC)
! Do k=1,Nz
!     j=Ny ! check if we applied on ghost layer or on the first fluid layer
!     Do i=1,Nx
!         l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
!         !top plane
!         f1(l,:)=2.d0*f1(l-Nxtotal,:)-f1(l-2*Nxtotal,:) !Check range is whole plane
!         f2(l,:)=2.d0*f2(l-Nxtotal,:)-f2(l-2*Nxtotal,:) !or only the gap-line for this
!         f5(l,:)=2.d0*f5(l-Nxtotal,:)-f5(l-2*Nxtotal,:) !extrapolation
!         f6(l,:)=2.d0*f6(l-Nxtotal,:)-f6(l-2*Nxtotal,:)
!         f4(l,:)=f1(l,:)
!         f3(l,:)=f2(l,:)
!         f8(l,:)=f5(l,:)
!         f7(l,:)=f6(l,:)
!     End do
! End do
! !$OMP END DO 

! !$OMP DO SCHEDULE(STATIC) 
! Do k=1,Nz
!     j=0 ! check if we applied on ghost layer or on the first fluid layer
!     Do i=1,Nx
!         l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
!         !bottom plane
!         f2(l,:)=f3(l+2*Nxtotal,:)
!         f1(l,:)=f4(l+2*Nxtotal,:)
!         f6(l,:)=f7(l+2*Nxtotal,:)
!         f5(l,:)=f8(l+2*Nxtotal,:)
!     End do
! End do
! !$OMP END DO NOWAIT
! 
! !$OMP DO SCHEDULE(STATIC)
! Do k=1,Nz
!     j=Ny+1 ! check if we applied on ghost layer or on the first fluid layer
!     Do i=1,Nx
!         l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
!         !top plane
!         f4(l,:)=f1(l-2*Nxtotal,:)
!         f3(l,:)=f2(l-2*Nxtotal,:)
!         f8(l,:)=f5(l-2*Nxtotal,:)
!         f7(l,:)=f6(l-2*Nxtotal,:)
!     End do
! End do
! !$OMP END DO 

!$OMP DO SCHEDULE(STATIC) 
Do k=1,Nz
    j=1 ! check if we applied on ghost layer or on the first fluid layer
    Do i=1,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        !bottom plane
        f2(l,:)=f3(l,:)
        f1(l,:)=f4(l,:)
        f6(l,:)=f7(l,:)
        f5(l,:)=f8(l,:)
    End do
End do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
Do k=1,Nz
    j=Ny ! check if we applied on ghost layer or on the first fluid layer
    Do i=1,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        !top plane
        f4(l,:)=f1(l,:)
        f3(l,:)=f2(l,:)
        f8(l,:)=f5(l,:)
        f7(l,:)=f6(l,:)
    End do
End do
!$OMP END DO 
!=======================================================================
!     Boundary condition on symmetric planes z=
!=======================================================================
! !$OMP DO SCHEDULE(STATIC)
! Do j=1,Ny
!     k=0 ! check if we applied on ghost layer or on the first fluid layer
!     Do i=1,Nx
!         l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
!         !back plane
!         f1(l,:)=f5(l+2*Nxytotal,:)
!         f2(l,:)=f6(l+2*Nxytotal,:)
!         f3(l,:)=f7(l+2*Nxytotal,:)
!         f4(l,:)=f8(l+2*Nxytotal,:)
!     End do
! End do
! !$OMP END DO NOWAIT
! 
! !$OMP DO SCHEDULE(STATIC)
! Do j=1,Ny
!     k=Nz+1 ! check if we applied on ghost layer or on the first fluid layer
!     Do i=1,Nx
!         l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
!         !front plane
!         f5(l,:)=f1(l-2*Nxytotal,:)
!         f6(l,:)=f2(l-2*Nxytotal,:)
!         f7(l,:)=f3(l-2*Nxytotal,:)
!         f8(l,:)=f4(l-2*Nxytotal,:)
!     End do
! End do
! !$OMP END DO

!$OMP DO SCHEDULE(STATIC)
Do j=1,Ny
    k=1 ! check if we applied on ghost layer or on the first fluid layer
    Do i=1,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        !back plane
        f1(l,:)=f5(l,:)
        f2(l,:)=f6(l,:)
        f3(l,:)=f7(l,:)
        f4(l,:)=f8(l,:)
    End do
End do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
Do j=1,Ny
    k=Nz ! check if we applied on ghost layer or on the first fluid layer
    Do i=1,Nx
        l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        !front plane
        f5(l,:)=f1(l,:)
        f6(l,:)=f2(l,:)
        f7(l,:)=f3(l,:)
        f8(l,:)=f4(l,:)
    End do
End do
!$OMP END DO

!=======================================================================
!     Macroscopic parameter evaluation
!=======================================================================
!------------------------------------------------------------------------
!           Rho, Ux, Uy
!------------------------------------------------------------------------
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
!------------------------------------------------------------------------
!           Flow rate
!------------------------------------------------------------------------
    iteration1=iteration1+1
    If (mod(iteration1,interval)==0) then
        mass2=0.d0
        Do k=2,Nz-1
        Do j=2,Ny-1
            l=(column+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
            mass2=mass2+Ux(l)*ds*ds
        Enddo
        Enddo
        k=1
        Do j=2,Ny-1
            l=(column+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
            mass2=mass2+(Ux(l)+Ux(l+(Nz-1)*Nxytotal))*ds*ds*0.5d0
        Enddo
        mass2=(mass2+0.25d0*ds*ds*(Ux(column+ghostLayer*Nxtotal+ghostLayer*Nxytotal) &
        &      +Ux(column+(ghostLayer+Ny-1)*Nxtotal+ghostLayer*Nxytotal) &
        &      +Ux(column+ghostLayer*Nxtotal+(ghostLayer+Nz-1)*Nxytotal) &
        &      +Ux(column+(ghostLayer+Ny-1)*Nxtotal+(ghostLayer+Nz-1)*Nxytotal) &
        & ))*sqrt(1.d0/2.d0)*4.d0/PressDrop
        error=dabs(1.d0-mass2/mass)/(interval)
        mass=mass2
        permeability=mass*seriesKn(iKn)*sqrt(4.d0/pi)

!------------------------------------------------------------------------
!           At the wall nodes
!------------------------------------------------------------------------
        write(*,"( 1I10, 3E14.4)")  iteration1,  mass,  permeability, error
    endif
Enddo

call cpu_time(end_time)
!***********************************************************************
!   Print out results
!***********************************************************************
    open(20,file='Field.dat',STATUS="REPLACE")
        write(20,*) ' TITLE=" Field"'
        write(20,*) ' VARIABLES=x,y,z,Rho,Ux,Uy,Uz'
        write(20,*)'ZONE T=final, I=',Nx,', J=',Ny,', K=',Nz,', F=POINT'
        Do k=1,Nz
        Do j=1,Ny
            Do i=1,Nx
                l=(i+ghostLayer)+(j+ghostLayer-1)*Nxtotal+(k+ghostLayer-1)*Nxytotal
        If (image(l)==fluid) then
                    write(20,'(5ES15.3E3)') (i-1)*ds, (j-1)*ds, (k-1)*ds, Rho(l)+1.d0, Ux(l), Uy(l), Uz(l)
                else
            write(20,'(5ES15.3E3)') (i-1)*ds, (j-1)*ds, (k-1)*ds, 0.d0, 0.d0, 0.d0, 0.d0
        endif
        Enddo
        Enddo
        Enddo
    close(20)

    open(22,file='Results.dat', position="append")
        write(22,"(4E14.4, 1I15, E15.3)") seriesKn(iKn), mass, permeability,error,iteration1,(end_time-start_time)/3600.d0
    close(22)

Enddo

DEALLOCATE(vecWall)
DEALLOCATE(dir1, dir2, dir3, dir4, dir5, dir6, dir7, dir8)
DEALLOCATE(coef1, coef2, coef3, coef4, coef5, coef6, coef7, coef8)
DEALLOCATE(f1,f2,f3,f4,f5,f6,f7,f8)
DEALLOCATE(f1w,f2w,f3w,f4w,f5w,f6w,f7w,f8w)
DEALLOCATE(Rho, Ux, Uy, Uz)
DEALLOCATE(image,which_corner)

END PROGRAM linearised_BGK_fix_point
