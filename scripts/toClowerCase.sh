#!/bin/bash

## This scripts is used for clean source code
## It will use the owerCase.txt .


for i in $(cat lowerCase.txt); do
    goodi=$(echo ${i}|tr '[:upper:]' '[:lower:]')
    #echo ${i} ${goodi}
    sed -i -e "s/${i}/${goodi}/g" ../src/*.F90
    sed -i -e "s/SHifT/shift/g" ../src/*.F90
    sed -i -e "s/doUBLE/DOUBLE/g" ../src/*.F90
    sed -i -e "s/COor/COOR/g" ../src/*.F90
    sed -i -e "s/WorLD/WORLD/g" ../src/*.F90
    sed -i -e "s/PARAfile/PARAFILE/g" ../src/*.F90
    sed -i -e "s/IMAGEfile/IMAGEFILE/g" ../src/*.F90
    sed -i -e "s/ERRor/ERROR/g" ../src/*.F90
    sed -i -e "s/MPI_status_SIZE/MPI_STATUS_SIZE/g" ../src/*.F90
    sed -i -e "s/ISend/ISEND/g" ../src/*.F90
    sed -i -e "s/OMP do/OMP DO/g" ../src/*.F90
    sed -i -e "s/OMP end/OMP END/g" ../src/*.F90
    sed -i -e "s/OMP enddo/OMP ENDDO/g" ../src/*.F90
    sed -i -e "s/OMP end do/OMP END DO/g" ../src/*.F90
done
