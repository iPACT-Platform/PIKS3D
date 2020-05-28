!/------------------------------------------------------------\
!                                                             |
!  $$$$$$$\ $$$$$$\ $$\   $$\  $$$$$$\   $$$$$$\  $$$$$$$\    |
!  $$  __$$\\_$$  _|$$ | $$  |$$  __$$\ $$ ___$$\ $$  __$$\   |
!  $$ |  $$ | $$ |  $$ |$$  / $$ /  \__|\_/   $$ |$$ |  $$ |  |
!  $$$$$$$  | $$ |  $$$$$  /  \$$$$$$\    $$$$$ / $$ |  $$ |  |
!  $$  ____/  $$ |  $$  $$<    \____$$\   \___$$\ $$ |  $$ |  |
!  $$ |       $$ |  $$ |\$$\  $$\   $$ |$$\   $$ |$$ |  $$ |  |
!  $$ |     $$$$$$\ $$ | \$$\ \$$$$$$  |\$$$$$$  |$$$$$$$  |  |
!  \__|     \______|\__|  \__| \______/  \______/ \_______/   |
!                                                             |
!       Parallel Image-based Kinetic Solver (3D)              |
!\------------------------------------------------------------/

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
! module    : velocityGrid
!-------------------------------------------------------------------------------
! This is a module for velocity space configurations of 3D DVM parallel solver. 
! For details:
!
! [1]   M.T. Ho, L. Zhu, L. Wu, P. Wang, Z. Guo, Z.-H. Li, Y. Zhang
!       "A multi-level parallel solver for rarefied gas flows in porous media"
! 		Computer Physics Communications, 234 (2019), pp. 14-25
!
!	See Section 2.2 of Ref.[1]
!-------------------------------------------------------------------------------

module velocityGrid
use gaussHermite

implicit none
save

! 2*Nc_fundamental is the number of the discrete velocities in one axis,
! or the order of the quadrature.
! The total number of discrete velocities is Nv=(2*Nc_fundamental)^3
! Nc_fundamental to be read from NML: velocityNml
integer :: Nc_fundamental

! Whether to use the half-range variant of Gauss-Hermite quadrature
logical :: halfRange

!Number of moleculer velocities in 3D space and a octant of it.
integer :: Nc, Nc8

!abscissae and weighting Hermite quadrature, dimension(Nc_fundamental)
double precision, dimension (:), allocatable :: xi, weight1D 

!molecular velocity and weighting, dimension(Nc)
double precision, dimension (:), allocatable :: cx, cy, cz, w

!specular wall's normal vector in X, Y direction, dimension(Nc)
integer, dimension (:), allocatable :: oppositeX, oppositeY, oppositeZ

!constant PI
double precision, parameter :: PI=datan(1.d0)*4.d0

!half range flux of discrete velocity grid 
double precision :: DiffFlux

contains
    !> @brief setup the discrete velocity grid
    !!
    !! Use Nc_fundamental[integer] and halfRange[logical],
    !! to construct the discrete velocity set and weight coefficients.
    !! If halfRange = true, then use half-range Gauss-Hermite quadrature
    !! points as discrete velocities, otherwise use the full-range Gauss-Hermite quadrature.
    subroutine setupVelocityGrid
        implicit none
        ! index used later
        integer :: l, m, n, k

        ! Nc_fundamental has been initialized from the NML input, so we know 
        ! the number of the discrete velocities in half of one axis
        allocate(xi(Nc_fundamental))
        allocate(weight1D(Nc_fundamental))

        ! Use Gauss-Hermite, assign the discrete velocities point and weights,
        ! Only the 2nd, 4th, 6th, 8th ... orders are supported
        if(.not. halfRange) then
            select case (Nc_fundamental)
                case(2)
                    xi = xi2
                    weight1D = wi2
                case(4)
                    xi = xi4
                    weight1D = wi4
                case(6)
                    xi = xi6
                    weight1D = wi6
                case(8)
                    xi = xi8
                    weight1D = wi8
                case(12)
                    xi = xi12
                    weight1D = wi12
                case(16)
                    xi = xi16
                    weight1D = wi16
                case default
                    print*, "Nc_fundamental xi/wi for ", Nc_fundamental, &
                     "has not been provided"
                    deallocate(xi)
                    deallocate(weight1D)
            endselect
        ! Use half-range Gauss-Hermite
        else
            select case (Nc_fundamental)
                case(2)
                    xi = hxi2
                    weight1D = hwi2
                case(4)
                    xi = hxi4
                    weight1D = hwi4
                case(6)
                    xi = hxi6
                    weight1D = hwi6
                case(8)
                    xi = hxi8
                    weight1D = hwi8
                case(12)
                    xi = hxi12
                    weight1D = hwi12
                case(16)
                    xi = hxi16
                    weight1D = hwi16
                case default
                    print*, "Nc_fundamental xi/wi for ", Nc_fundamental, &
                     "has not been provided"
                    deallocate(xi)
                    deallocate(weight1D)
            endselect
        endif

        ! Number of the discrete velocity in 3-dimesions
        Nc=(2*Nc_fundamental)**3
        ! Number of dimcrete velocities in one quadrants
        Nc8 = Nc/8

        ! The discrete velocities components
        allocate(cx(Nc8), cy(Nc8), cz(Nc8), w(Nc8))
        allocate(oppositeX(Nc), oppositeY(Nc), oppositeZ(Nc))

        ! BUG find 2017-06-09
        ! should be properly scaled
        do l=1,Nc_fundamental
            xi(l) = xi(l)*dsqrt(2.d0)
        enddo
        
        ! Set the components of the discrete velocities and the total weights
        do l=1,Nc_fundamental
            do m=1,Nc_fundamental
                do n=1,Nc_fundamental
                    k=n+(m-1)*Nc_fundamental+(l-1)*Nc_fundamental**2
                    cx(k) = xi(n)
                    cy(k) = xi(m)
                    cz(k) = xi(l)
                    w(k) = weight1D(n)*weight1D(m)*weight1D(l)
                enddo
            enddo
        enddo  
        
        ! diffFlux is a velocity-set related constant, see denominator of Eq.(9) in Ref.[1]
        diffFlux=0.d0
        do l=1,Nc8
            DiffFlux=DiffFlux+cz(l)*w(l)
        enddo
        DiffFlux=DiffFlux*4

    end subroutine setupVelocityGrid

end module velocityGrid
