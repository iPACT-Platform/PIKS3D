program main
implicit none

integer :: Nxtotal = 38 
integer :: Nytotal = 34
integer :: Nztotal = 34
integer :: k = 20
integer :: j = 19
integer :: i = 24
integer :: xlg, xug, ylg, yug, zlg, zug

integer :: l, kk, jj, ii
integer :: Nxytotal
Nxytotal = Nxtotal*Nytotal
xlg = -1
ylg = -1
zlg = -1

l = (k-zlg)*Nxtotal*Nytotal + (j-ylg)*Nxtotal + (i-xlg+1)

!kk = mod(l,Nxytotal) + zlg ! to be checked
!jj = mod(l-(kk-zlg)*Nxytotal, Nxtotal) + ylg
!ii = l-(kk-zlg)*Nxytotal-(jj-ylg)*Nxtotal + xlg -1 
kk = l/Nxytotal + zlg ! to be checked
jj = (l-(kk-zlg)*Nxytotal)/Nxtotal + ylg
ii = l-(kk-zlg)*Nxytotal-(jj-ylg)*Nxtotal + xlg -1 
print*, "l=", l
print*, kk, jj, ii
end program main