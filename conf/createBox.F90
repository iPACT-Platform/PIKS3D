program main
implicit none
integer, parameter :: NX = 50
integer, parameter :: NY = 50
integer, parameter :: NZ = 50
integer, parameter :: startI = 12
integer, parameter :: endI   = 39
character, allocatable, dimension (:,:,:) :: flag

integer :: i, j, k
integer :: openstatus
double precision :: disSq

allocate(flag(NX, NY, NZ))

do k = 1, NZ
    do j = 1, NY
        do i = 1, NX
            if ( (i > startI .and. i < endI) &
                 (j > startI .and. j < endI) &
                 (k > startI .and. k < endI) &) then
                flag(i,j,k) = char(1)
            else
                flag(i,j,k) = char(0)
            endif
        enddo
    enddo
enddo

open(unit=11, file = 'singleBox.raw', status='replace', &
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
