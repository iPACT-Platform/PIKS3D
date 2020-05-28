!-----------------------------------------------------------
! This program creates an simple binary geometry flag file 
! which can used as the `rawImageFile` for the program PIKS3D.
! It create a solid sphere in a box of dimension 100^3 such 
! that the sphere locates on the center of the box, and the 
! void space takes 75% of the box's total volume.
! From this program, the users of PIKS3D will find out the 
! format of the binary image file definging a geometry
!-----------------------------------------------------------

program main
implicit none
integer, parameter :: NX = 100
integer, parameter :: NY = 100
integer, parameter :: NZ = 100
double precision, parameter :: porosity = 0.75d0
integer, parameter :: obstR = floor((NX-1)*((1.d0-porosity)*3.d0/4.d0/(datan(1.d0)*4.d0))**(1.d0/3.d0))
double precision, parameter :: obstX = (NX+1)/2.d0
double precision, parameter :: obstY = 1.d0
double precision, parameter :: obstZ = 1.d0
character, allocatable, dimension (:,:,:) :: flag


integer :: i, j, k
integer :: openstatus
double precision :: disSq

allocate(flag(NX, NY, NZ))

do k = 1, NZ
    do j = 1, NY
        do i = 1, NX
            disSq = (k - obstZ)**2 + (j - obstY)**2 + (i - obstX)**2
            if (disSq < obstR**2) then
                flag(i,j,k) = char(1)
            else
                flag(i,j,k) = char(0)
            endif
        enddo
    enddo
enddo

open(unit=11, file = 'singleSphere.raw', status='replace', &
    form='unformatted', action='write', access='stream', &
    iostat=openstatus)

if(openstatus /= 0) then
    print *, 'Could not open flag file'
    stop
endif
write(11) flag
close(11)
print *, 'write file done'
deallocate(flag)
end program
