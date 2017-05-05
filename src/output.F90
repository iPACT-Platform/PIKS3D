module output
use flow
use physicalGrid
use mpiParams
use solver

implicit none
save

contains
    subroutine chkConverge
        implicit none
        include "mpif.h"

        ! local vars
        integer :: byl, byu, bzl, bzu, k, j, i, l, MPI_ERR
        double precision :: massInner, massNoth, massSuth, massFrnt, massBack
        double precision :: massNF, massNB, massSF, massSB
        double precision :: massLocal, mass2

        massInner = 0.d0
        massNoth = 0.d0
        massSuth = 0.d0
        massFrnt = 0.d0
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
        if(yu == ymax) byu = yu - 1 !if most north block
        if(zl == zmin) bzl = zl + 1 !if most back block
        if(zu == zmax) bzu = zu - 1 !if most frnt block

        !mass2=0.d0
        if(xl == xmin) then !only left most processors
            do k=bzl,bzu
                do j=byl,byu
                    l=(k-zlg)*Nxytotal + (j-ylg)*Nxtotal + column+ghostLayers
                    massInner=massInner+Ux(l)*ds*ds
                enddo
            enddo
            if (yl == ymin) then !only south most processors
                do k=bzl,bzu
                    l = (k-zlg)*Nxytotal + ghostLayers*Nxtotal + column+ghostLayers
                    massSuth = massSuth + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (yu == ymax) then !only north most processors
                do k=bzl,bzu
                    l = (k-zlg)*Nxytotal + (Nysub+ghostLayers-1)*Nxtotal + column+ghostLayers
                    massNoth = massNoth + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (zl == zmin) then !only back most processors
                do j=byl,byu
                    l = ghostLayers*Nxytotal + (j-ylg)*Nxtotal + column+ghostLayers
                    massBack = massBack + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (zu == zmax) then !only frnt most processors
                do j=byl,byu
                    l = (Nzsub+ghostLayers-1)*Nxytotal + (j-ylg)*Nxtotal + column+ghostLayers
                    massFrnt = massFrnt + Ux(l)*ds*ds*0.5d0
                enddo
            endif
            if (yl==ymin.and.zl==zmin) then
                l = ghostLayers*Nxytotal + ghostLayers*Nxtotal + column+ghostLayers
                massSB = Ux(l)*ds*ds*0.25d0
            endif
            if (yl==ymin.and.zu==zmax) then
                l = (Nzsub+ghostLayers-1)*Nxytotal + ghostLayers*Nxtotal + column+ghostLayers
                massSF = Ux(l)*ds*ds*0.25d0
            endif
            if (yu==ymax.and.zl==zmin) then
                l = ghostLayers*Nxytotal + (Nysub+ghostLayers-1)*Nxtotal + column+ghostLayers
                massNB = Ux(l)*ds*ds*0.25d0
            endif
            if (yu==ymax.and.zu==zmax) then
                l = (Nzsub+ghostLayers-1)*Nxytotal + (Nysub+ghostLayers-1)*Nxtotal + column+ghostLayers
                massNF = Ux(l)*ds*ds*0.25d0
            endif


            ! debug
            massLocal = (massInner + massSuth + massNoth + massFrnt + massBack &
                + massSB + massSF + massNB + massNF) &
                /PressDrop/(1.d0/Ref_L)**2
                !+ massSB + massSF + massNB + massNF)*dsqrt(1.d0/2.d0)*4.d0/1.d0

            ! reduction
            call MPI_ALLREDUCE(massLocal, mass2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               mpi_comm_inlet, MPI_ERR)

            !PRINT*, "mass = ", mass
            !DEBUG
            error=dabs(1.d0-mass2/mass)/(chkConvergeStep)

            mass=mass2
            if (proc == master) then           
                permeability=mass*Kn*sqrt(4.d0/PI)  
                write(*,"( 1I10, 3ES15.6)")  iStep,  mass,  permeability, error
                ! open(22,file='Results.dat', position="append")
                ! write(22,'(4ES15.6, 1I15)') Kn, mass, permeability, error, iStep
                ! close(22)
            endif
        endif ! xl==xmin
        !bcast error so every process in WORLD can stop
        CALL MPI_BCAST(error, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_VGRID, MPI_ERR)
    end subroutine chkConverge

    subroutine saveFlowField
        integer :: j, i, k, l
        character(13) fname
        character(9) ztitle
        write(fname, '(A, I0.3, A)') 'Field_', proc, '.dat'
        write(ztitle, '(A, I0.3,A)') '"proc', proc, '"'
        open(20,file=fname,STATUS="REPLACE")
        write(20,*) ' TITLE=" Field"'
        write(20,*) ' VARIABLES=x,y,z,flag,Rho,Ux,Uy,Uz'
        write(20,'(A,A,A,I,A,I,A,I,A)') ' ZONE T=', ztitle, ', I=', Nxsub,', J=', Nysub, ', K=', Nzsub, ' F=POINT'
        Do k=zl,zu
            Do j=yl,yu
                Do i=xl,xu
                    l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                    If (image(i,j,k)==fluid) then
                        write(20,'(4I10,4ES15.6)') i, j, k, 0, Rho(l)+1.d0, Ux(l), Uy(l), Uz(l)
                    else
                        write(20,'(4I10,4ES15.6)') i, j, k, 1, 0.d0, 0.d0, 0.d0, 0.d0
                    Endif   
                Enddo
            Enddo
        Enddo
        close(20)
    end subroutine saveFlowField


    SUBROUTINE saveFlowFieldVTK
        IMPLICIT NONE
        INTEGER :: i, j, k, l,MPI_ERR, IO_ERR
        character(13) fname

        write(fname, '(A, I0.3, A)') 'Field_', proc, '.vtk'
        OPEN(UNIT = 12, FILE = fname, STATUS = "REPLACE", POSITION = "APPEND", &
          IOSTAT = IO_ERR)
        IF ( IO_ERR == 0 ) THEN
            WRITE(12,'(A)')"# vtk DataFile Version 2.0"
            WRITE(12,'(A)')"DVM MPI"
            WRITE(12,'(A)')"ASCII"
            WRITE(12,*)
            WRITE(12,'(A)')"DATASET STRUCTURED_POINTS"
            WRITE(12,*)"DIMENSIONS",Nxsub,Nysub,Nzsub
            WRITE(12,*)"ORIGIN",xl,yl,zl
            WRITE(12,*)"SPACING",1,1,1
            WRITE(12,*)
            WRITE(12,*)"POINT_DATA",Nxsub*Nysub*Nzsub
            WRITE(12,*)
            WRITE(12,'(A)')"SCALARS Rho double"
            WRITE(12,'(A)')"LOOKUP_TABLE default"
            DO k = zl, zu
                DO j = yl, yu
                    DO i = xl, xu
                        l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                        If (image(i,j,k)==fluid) then
                            write(12,'(ES15.6)') Rho(l)+1.d0
                        else
                            write(12,'(ES15.6)') 0.d0
                        Endif   
                    END DO
                END DO
            END DO
            WRITE(12,*)
            WRITE(12,'(A)')"VECTORS Velocity double"
            DO k = zl, zu
               DO j = yl, yu
                  DO i = xl, xu
                     l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                     If (image(i,j,k)==fluid) then
                         write(12,'(3ES15.6)') Ux(l), Uy(l), Uz(l)
                     else
                         write(12,'(3ES15.6)') 0.d0, 0.d0, 0.d0
                     Endif  
                  END DO
               END DO
            END DO
            CLOSE(UNIT = 12)
        ELSE
            CALL memFree
            CALL MPI_FINALIZE(MPI_ERR)
            STOP "Error: Unable to open output vtk file."
        END IF

        RETURN
    END SUBROUTINE saveFlowFieldVTK

    SUBROUTINE saveFlowFieldVTI
        INTEGER :: i, j, k, l,MPI_ERR, IO_ERR
        character(13) fname
        character(10) pfname
        integer :: exl, exu, eyl, eyu, ezl, ezu
        integer :: wexl, wexu, weyl, weyu, wezl, wezu
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

        write(fname, '(A, I0.3, A)') 'Field_', proc, '.vti'
        OPEN(UNIT = 13, FILE = fname, STATUS = "REPLACE", POSITION = "APPEND", &
          IOSTAT = IO_ERR)
        IF ( IO_ERR == 0 ) THEN
            WRITE(13,'(A)') '<?xml version="1.0"?>'
            WRITE(13,'(A)') '<VTKFile type="ImageData">'
            WRITE(13,'(A, 6I4, A)') '<ImageData WholeExtent="', exl, exu, & 
                eyl, eyu, ezl, ezu, ' " Origin="0 0 0" Spacing="1 1 1">'
            WRITE(13, '(A, 6I4, A)') '<Piece Extent="', exl, exu, eyl, eyu, &
                ezl, ezu, '">'
            WRITE(13, '(A)') '<CellData Scalars="flag Rho" Vectors ="U">'
            WRITE(13, '(A)') '<DataArray type="Int32" Name="flag" format="ascii">'
            DO k = zl, zu
                DO j = yl, yu
                    DO i = xl, xu
                        write(13,'(I)') image(i,j,k)
                    END DO
                END DO
            END DO
            WRITE(13, '(A)') '</DataArray>'
            WRITE(13, '(A)') '<DataArray type="Float32" Name="Rho" format="ascii">'
            DO k = zl, zu
                DO j = yl, yu
                    DO i = xl, xu
                        l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                        If (image(i,j,k)==fluid) then
                            write(13,'(ES15.6)') Rho(l)+1.d0
                        else
                            write(13,'(ES15.6)') 0.d0
                        Endif   
                    END DO
                END DO
            END DO
            WRITE(13, '(A)') '</DataArray>'
            WRITE(13, '(A)') '<DataArray type="Float32" Name="U" format="ascii" &
                & NumberOfComponents="3">'
            DO k = zl, zu
               DO j = yl, yu
                  DO i = xl, xu
                     l= (k-zlg)*Nxytotal + (j-ylg)*Nxtotal + i-xlg+1
                     If (image(i,j,k)==fluid) then
                         write(13,'(3ES15.6)') Ux(l), Uy(l), Uz(l)
                     else
                         write(13,'(3ES15.6)') 0.d0, 0.d0, 0.d0
                     Endif  
                  END DO
               END DO
            END DO
            WRITE(13, '(A)') '</DataArray>'
            WRITE(13, '(A)') '</CellData>'
            WRITE(13, '(A)') '</Piece>'
            WRITE(13, '(A)') '</ImageData>'
            WRITE(13, '(A)') '</VTKFile>'
            CLOSE(UNIT = 13)
        ELSE
            CALL memFree
            CALL MPI_FINALIZE(MPI_ERR)
            STOP "Error: Unable to open output vti file."
        END IF


        if (proc == master) then  
            write(pfname, '(A)') 'Field.pvti'
            OPEN(UNIT = 14, FILE = pfname, STATUS = "REPLACE", POSITION = "APPEND", &
              IOSTAT = IO_ERR)
            WRITE(14,'(A)') '<?xml version="1.0"?>'
            WRITE(14,'(A)') '<VTKFile type="PImageData">'
            WRITE(14,'(A, 6I4, A)') '<PImageData WholeExtent="', wexl, wexu, & 
                weyl, weyu, wezl, wezu, ' " Origin="0 0 0" Spacing="1 1 1">'
            WRITE(14,'(A)') '<PCellData Scalars="flag Rho" Vectors ="U">'
            WRITE(14, '(A)') '<DataArray type="Int32" Name="flag" format="ascii"/>'
            WRITE(14, '(A)') '<DataArray type="Float32" Name="Rho" format="ascii"/>'
            WRITE(14, '(A)') '<DataArray type="Float32" Name="U" format="ascii" &
                & NumberOfComponents="3"/>'
            WRITE(14, '(A)') '</PCellData>'
            do l=1, nprocs
                exl = sub_ext(1,l) - 1
                exu = sub_ext(2,l)
                eyl = sub_ext(3,l) - 1
                eyu = sub_ext(4,l)
                ezl = sub_ext(5,l) - 1
                ezu = sub_ext(6,l)
                write(fname, '(A, I0.3, A)') 'Field_', l-1, '.vti'
                WRITE(14, '(A, 6I4, A)') '<Piece Extent="', exl, exu, eyl, eyu, &
                    ezl, ezu, '" Source="'//fname//'"/>'
            enddo
            WRITE(14, '(A)') '</PImageData>'
            WRITE(14, '(A)') '</VTKFile>'
            CLOSE(UNIT=14)
        endif

    END SUBROUTINE saveFlowFieldVTI


    subroutine memFree
        deallocate (array3D, array3Dg)
        deallocate (image, which_corner, vecWall)
        deallocate (dir1, dir2, dir3, dir4, dir5, dir6, dir7, dir8)
        deallocate (coef1, coef2, coef3, coef4, coef5, coef6, coef7, coef8)
        deallocate (f1,f2,f3,f4,f5,f6,f7,f8)
        !deallocate (f1w,f2w,f3w,f4w,f5w,f6w,f7w,f8w)
        deallocate (Rho, Ux, Uy, Uz)
        deallocate (f_west_snd, f_west_rcv, f_east_snd, f_east_rcv)
        deallocate (f_noth_snd, f_noth_rcv, f_suth_snd, f_suth_rcv)
        deallocate (f_back_snd, f_back_rcv, f_frnt_snd, f_frnt_rcv)
    end subroutine memFree
end module output
