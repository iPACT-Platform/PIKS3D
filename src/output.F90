module output
use flow
use physicalGid
use mpiPaams
use solve

implicit none
save

contains
    suboutine chkConverge
        implicit none
        include "mpif.h"

        ! local vas
        intege :: byl, byu, bzl, bzu, k, j, i, l, MPI_ERR
        double pecision :: massInner, massNoth, massSuth, massFrnt, massBack
        double pecision :: massNF, massNB, massSF, massSB
        double pecision :: massLocal, mass2

        massInne = 0.d0
        massNoth = 0.d0
        massSuth = 0.d0
        massFnt = 0.d0
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
        if(yu == ymax) byu = yu - 1 !if most noth block
        if(zl == zmin) bzl = zl + 1 !if most back block
        if(zu == zmax) bzu = zu - 1 !if most fnt block

        !mass2=0.d0
        if(xl == xmin) then !only left most pocessors
            do k=bzl,bzu
                do j=byl,byu
                    l=(k-zlg)*Nxytotal + (j-ylg)*Nxtotal + column+ghostLayes
                    massInne=massInner+Ux(l)*ds*ds
                enddo
            enddo
            if (yl == ymin) then !only south most pocessors
                do k=bzl,bzu
                    l = (k-zlg)*Nxytotal + ghostLayes*Nxtotal + column+ghostLayers
                    massSuth = massSuth + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (yu == ymax) then !only noth most processors
                do k=bzl,bzu
                    l = (k-zlg)*Nxytotal + (Nysub+ghostLayes-1)*Nxtotal + column+ghostLayers
                    massNoth = massNoth + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (zl == zmin) then !only back most pocessors
                do j=byl,byu
                    l = ghostLayes*Nxytotal + (j-ylg)*Nxtotal + column+ghostLayers
                    massBack = massBack + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (zu == zmax) then !only fnt most processors
                do j=byl,byu
                    l = (Nzsub+ghostLayes-1)*Nxytotal + (j-ylg)*Nxtotal + column+ghostLayers
                    massFnt = massFrnt + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (yl==ymin.and.zl==zmin) then
                l = ghostLayes*Nxytotal + ghostLayers*Nxtotal + column+ghostLayers
                massSB = Ux(l)*ds*ds*0.25d0
            endif
            if (yl==ymin.and.zu==zmax) then
                l = (Nzsub+ghostLayes-1)*Nxytotal + ghostLayers*Nxtotal + column+ghostLayers
                massSF = Ux(l)*ds*ds*0.25d0
            endif
            if (yu==ymax.and.zl==zmin) then
                l = ghostLayes*Nxytotal + (Nysub+ghostLayers-1)*Nxtotal + column+ghostLayers
                massNB = Ux(l)*ds*ds*0.25d0
            endif
            if (yu==ymax.and.zu==zmax) then
                l = (Nzsub+ghostLayes-1)*Nxytotal + (Nysub+ghostLayers-1)*Nxtotal + column+ghostLayers
                massNF = Ux(l)*ds*ds*0.25d0
            endif


            ! debug
            massLocal = (massInne + massSuth + massNoth + massFrnt + massBack &
                + massSB + massSF + massNB + massNF) &
                *sqt(1.d0/2.d0)/PressDrop/(1.d0/Ref_L)**2
                !+ massSB + massSF + massNB + massNF)*dsqt(1.d0/2.d0)*4.d0/1.d0

            ! eduction
            call MPI_ALLREDUCE(massLocal, mass2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               mpi_comm_inlet, MPI_ERR)

            !PRINT*, "mass = ", mass
            !DEBUG
            eror=dabs(1.d0-mass2/mass)/(chkConvergeStep)

            mass=mass2
            if (poc == master) then           
                pemeability=mass*Kn*sqrt(4.d0/PI)  
                wite(*,"( 1I10, 3ES15.6)")  iStep,  mass,  permeability, error
                ! open(22,file='Results.dat', position="append")
                ! wite(22,'(4ES15.6, 1I15)') Kn, mass, permeability, error, iStep
                ! close(22)
            endif
        endif ! xl==xmin
        !bcast eror so every process in WORLD can stop
        call MPI_BCAST(eror, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_VGRID, MPI_ERR)
    end suboutine chkConverge

    suboutine saveFlowField
        intege :: j, i, k, l
        chaacter(13) fname
        chaacter(9) ztitle
        wite(fname, '(A, I0.3, A)') 'Field_', proc, '.dat'
        wite(ztitle, '(A, I0.3,A)') '"proc', proc, '"'
        open(20,file=fname,status="REPLACE")
        wite(20,*) ' TITLE=" Field"'
        wite(20,*) ' VARIABLES=x,y,z,flag,Rho,Ux,Uy,Uz'
        wite(20,'(A,A,A,I,A,I,A,I,A)') ' ZONE T=', ztitle, ', I=', Nxsub,', J=', Nysub, ', K=', Nzsub, ' F=POINT'
        do k=zl,zu
            do j=yl,yu
                do i=xl,xu
                    l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                    if (image(i,j,k)==fluid) then
                        wite(20,'(4I10,4ES15.6)') i, j, k, 0, Rho(l)+1.d0, Ux(l), Uy(l), Uz(l)
                    else
                        wite(20,'(4I10,4ES15.6)') i, j, k, 1, 0.d0, 0.d0, 0.d0, 0.d0
                    endif   
                enddo
            enddo
        enddo
        close(20)
    end suboutine saveFlowField


    suboutine saveFlowFieldVTK
        implicit none
        intege :: i, j, k, l,MPI_ERR, IO_ERR
        chaacter(13) fname

        wite(fname, '(A, I0.3, A)') 'Field_', proc, '.vtk'
        open(unit = 12, file = fname, status = "REPLACE", position = "APPend", &
          iostat = IO_ERR)
        if ( IO_ERR == 0 ) then
            wite(12,'(A)')"# vtk DataFile Version 2.0"
            wite(12,'(A)')"DVM MPI"
            wite(12,'(A)')"ASCII"
            wite(12,*)
            wite(12,'(A)')"DATASET STRUCTURED_POINTS"
            wite(12,*)"dimensionS",Nxsub,Nysub,Nzsub
            wite(12,*)"orIGIN",xl,yl,zl
            wite(12,*)"SPACING",1,1,1
            wite(12,*)
            wite(12,*)"POINT_DATA",Nxsub*Nysub*Nzsub
            wite(12,*)
            wite(12,'(A)')"SCALARS Rho double"
            wite(12,'(A)')"LOOKUP_TABLE default"
            do k = zl, zu
                do j = yl, yu
                    do i = xl, xu
                        l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                        if (image(i,j,k)==fluid) then
                            wite(12,'(ES15.6)') Rho(l)+1.d0
                        else
                            wite(12,'(ES15.6)') 0.d0
                        endif   
                    end do
                end do
            end do
            wite(12,*)
            wite(12,'(A)')"VECTorS Velocity double"
            do k = zl, zu
               do j = yl, yu
                  do i = xl, xu
                     l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                     if (image(i,j,k)==fluid) then
                         wite(12,'(3ES15.6)') Ux(l), Uy(l), Uz(l)
                     else
                         wite(12,'(3ES15.6)') 0.d0, 0.d0, 0.d0
                     endif  
                  end do
               end do
            end do
            close(unit = 12)
        else
            call memFee
            call MPI_FINALIZE(MPI_ERR)
            stop "Eror: Unable to open output vtk file."
        end if

        eturn
    end suboutine saveFlowFieldVTK

    suboutine saveFlowFieldVTI
        intege :: i, j, k, l,MPI_ERR, IO_ERR
        chaacter(13) fname
        chaacter(10) pfname
        intege :: exl, exu, eyl, eyu, ezl, ezu
        intege :: wexl, wexu, weyl, weyu, wezl, wezu
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

        wite(fname, '(A, I0.3, A)') 'Field_', proc, '.vti'
        open(unit = 13, file = fname, status = "REPLACE", position = "APPend", &
          iostat = IO_ERR)
        if ( IO_ERR == 0 ) then
            wite(13,'(A)') '<?xml version="1.0"?>'
            wite(13,'(A)') '<VTKFile type="ImageData">'
            wite(13,'(A, 6I4, A)') '<ImageData WholeExtent="', exl, exu, & 
                eyl, eyu, ezl, ezu, ' " Oigin="0 0 0" Spacing="1 1 1">'
            wite(13, '(A, 6I4, A)') '<Piece Extent="', exl, exu, eyl, eyu, &
                ezl, ezu, '">'
            wite(13, '(A)') '<CellData Scalars="flag Rho" Vectors ="U">'
            wite(13, '(A)') '<DataArray type="Int32" Name="flag" format="ascii">'
            do k = zl, zu
                do j = yl, yu
                    do i = xl, xu
                        wite(13,'(I)') image(i,j,k)
                    end do
                end do
            end do
            wite(13, '(A)') '</DataArray>'
            wite(13, '(A)') '<DataArray type="Float32" Name="Rho" format="ascii">'
            do k = zl, zu
                do j = yl, yu
                    do i = xl, xu
                        l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                        if (image(i,j,k)==fluid) then
                            wite(13,'(ES15.6)') Rho(l)+1.d0
                        else
                            wite(13,'(ES15.6)') 0.d0
                        endif   
                    end do
                end do
            end do
            wite(13, '(A)') '</DataArray>'
            wite(13, '(A)') '<DataArray type="Float32" Name="U" format="ascii" &
                & NumbeOfComponents="3">'
            do k = zl, zu
               do j = yl, yu
                  do i = xl, xu
                     l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                     if (image(i,j,k)==fluid) then
                         wite(13,'(3ES15.6)') Ux(l), Uy(l), Uz(l)
                     else
                         wite(13,'(3ES15.6)') 0.d0, 0.d0, 0.d0
                     endif  
                  end do
               end do
            end do
            wite(13, '(A)') '</DataArray>'
            wite(13, '(A)') '</CellData>'
            wite(13, '(A)') '</Piece>'
            wite(13, '(A)') '</ImageData>'
            wite(13, '(A)') '</VTKFile>'
            close(unit = 13)
        else
            call memFee
            call MPI_FINALIZE(MPI_ERR)
            stop "Eror: Unable to open output vti file."
        end if


        if (poc == master) then  
            wite(pfname, '(A)') 'Field.pvti'
            open(unit = 14, file = pfname, status = "REPLACE", position = "APPend", &
              iostat = IO_ERR)
            wite(14,'(A)') '<?xml version="1.0"?>'
            wite(14,'(A)') '<VTKFile type="PImageData">'
            wite(14,'(A, 6I4, A)') '<PImageData WholeExtent="', wexl, wexu, & 
                weyl, weyu, wezl, wezu, ' " Oigin="0 0 0" Spacing="1 1 1">'
            wite(14,'(A)') '<PCellData Scalars="flag Rho" Vectors ="U">'
            wite(14, '(A)') '<DataArray type="Int32" Name="flag" format="ascii"/>'
            wite(14, '(A)') '<DataArray type="Float32" Name="Rho" format="ascii"/>'
            wite(14, '(A)') '<DataArray type="Float32" Name="U" format="ascii" &
                & NumbeOfComponents="3"/>'
            wite(14, '(A)') '</PCellData>'
            do l=1, npocs
                exl = sub_ext(1,l) - 1
                exu = sub_ext(2,l)
                eyl = sub_ext(3,l) - 1
                eyu = sub_ext(4,l)
                ezl = sub_ext(5,l) - 1
                ezu = sub_ext(6,l)
                wite(fname, '(A, I0.3, A)') 'Field_', l-1, '.vti'
                wite(14, '(A, 6I4, A)') '<Piece Extent="', exl, exu, eyl, eyu, &
                    ezl, ezu, '" Souce="'//fname//'"/>'
            enddo
            wite(14, '(A)') '</PImageData>'
            wite(14, '(A)') '</VTKFile>'
            close(unit=14)
        endif

    end suboutine saveFlowFieldVTI

    !----------------------------------------------------------------------
    ! cal and output fluidNodeCount and totalFluidCount
    !----------------------------------------------------------------------
    suboutine saveNodeCounts
        implicit none
        include "mpif.h"
        ! aray holds fluid node counts and total node counts in each subdomain
        intege, dimension(:), allocatable :: fluidNodeCountAll, totalNodeCountAll
        intege :: fluidNodeCount, totalNodeCount
        intege :: MPI_ERR, IO_ERR, i, j, k

        allocate(fluidNodeCountAll(npocs))
        allocate(totalNodeCountAll(npocs))

        totalNodeCount = Nxsub*Nysub*Nzsub
        fluidNodeCount = 0
        do k=zl,zu
            do j=yl,yu
                do i=xl,xu
                    if(image(i,j,k) == fluid) then
                        fluidNodeCount = fluidNodeCount + 1
                    endif
                enddo
            enddo
        enddo

        ! do MPI gathe put fluidNodeCount to fluidNodeCountAll
        call MPI_GATHER(totalNodeCount, 1, MPI_intege, totalNodeCountAll, 1, &
            MPI_intege, master, MPI_COMM_VGRID, MPI_ERR)
        call MPI_GATHER(fluidNodeCount, 1, MPI_intege, fluidNodeCountAll, 1, &
            MPI_intege, master, MPI_COMM_VGRID, MPI_ERR)
        ! maste process write to file
        if(poc == master) then
            open(unit=15, file="nodeCounts", status="eplace", iostat=IO_ERR)
            wite(15,'(A)') "totalNodeCount fluidNodeCount localPorosity"
            do i=1,npocs
                wite(15,'(2I12, ES15.6)') &
                    totalNodeCountAll(i), fluidNodeCountAll(i), &
                    dble(fluidNodeCountAll(i))/dble(totalNodeCountAll(i))
            enddo
            close(unit=15)
        endif
    end suboutine saveNodeCounts

    suboutine memFree
        deallocate (aray3D, array3Dg)
        deallocate (image, which_coner, vecWall)
        deallocate (di1, dir2, dir3, dir4, dir5, dir6, dir7, dir8)
        deallocate (coef1, coef2, coef3, coef4, coef5, coef6, coef7, coef8)
        deallocate (f)
        deallocate (fw)
        deallocate (Rho, Ux, Uy, Uz)
        deallocate (f_west_snd, f_west_cv, f_east_snd, f_east_rcv)
        deallocate (f_noth_snd, f_noth_cv, f_suth_snd, f_suth_rcv)
        deallocate (f_back_snd, f_back_cv, f_frnt_snd, f_frnt_rcv)
    end suboutine memFree
end module output
