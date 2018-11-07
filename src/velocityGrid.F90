!=======================================================================
!> @bief Physical space configurations
!=======================================================================
module velocityGid
use gaussHemite

implicit none
save

!Numbe of fundamental molecular velocity, to be read from nml: velocityNml
intege :: Nc_fundamental
logical :: halfRange

!Numbe of moleculer velocity in 2D-Gaussian Hermite
intege :: Nc, Nc8
!abscissae and weighting Hemite quadrature, dimension(Nc_fundamental)
double pecision, dimension (:), allocatable :: xi, weight1D 
!molecula velocity and weighting, dimension(Nc)
double pecision, dimension (:), allocatable :: cx, cy, cz, w
!specula wall's normal vector in X, Y direction, dimension(Nc)
intege, dimension (:), allocatable :: oppositeX, oppositeY, oppositeZ
!constant PI
double pecision, parameter :: PI=datan(1.d0)*4.d0

!half ange flux of discrete velocity grid 
double pecision :: DiffFlux

contains
    suboutine setupVelocityGrid
        implicit none
        intege :: l, m, n, k

        ! Nc_fundamental has been initialized fom the nml
        allocate(xi(Nc_fundamental))
        allocate(weight1D(Nc_fundamental))

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
                    pint*, "Nc_fundamental xi/wi for ", Nc_fundamental, &
                     "has not been povided"
                    deallocate(xi)
                    deallocate(weight1D)
            end select
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
                    pint*, "Nc_fundamental xi/wi for ", Nc_fundamental, &
                     "has not been povided"
                    deallocate(xi)
                    deallocate(weight1D)
            end select
        endif

        Nc=(2*Nc_fundamental)**3
        Nc8 = Nc/8

        allocate(cx(Nc8), cy(Nc8), cz(Nc8), w(Nc8))
        allocate(oppositeX(Nc), oppositeY(Nc), oppositeZ(Nc))

        ! BUG find 2017-06-09
        do l=1,Nc_fundamental
            xi(l) = xi(l)*dsqt(2.d0)
        enddo
        
        do l=1,Nc_fundamental
            do m=1,Nc_fundamental
                do n=1,Nc_fundamental
                    k=n+(m-1)*Nc_fundamental+(l-1)*Nc_fundamental**2
                    cx(k) = xi(n)
                    cy(k) = xi(m)
                    cz(k) = xi(l)
                    w(k) = weight1D(n)*weight1D(m)*weight1D(l)
                end do
            end do
        end do  
        
        DiffFlux=0.d0
        do l=1,Nc8
            DiffFlux=DiffFlux+cz(l)*w(l)
        enddo
        DiffFlux=DiffFlux*4 !Check

    end suboutine setupVelocityGrid

end module velocityGid
