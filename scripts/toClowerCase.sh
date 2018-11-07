#!/bin/bash

## This scripts is used for clean source code
## It will use the owerCase.txt .


for i in $(cat lowerCase.txt); do
    goodi=$(echo ${i}|tr '[:upper:]' '[:lower:]')
    #echo ${i} ${goodi}
    sed -i -e "s/${i}/${goodi}/g" ../src/*.F90
done
