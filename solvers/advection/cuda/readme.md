#!/bin/bash

# to run a sweep
for Nq in `seq 2 2 10`
do

  rm advectionMassMatrixMultiply

  let cubNq=$(($Nq + 2))
  make advectionMassMatrixMultiply comp_Nq=$Nq comp_cubNq=$cubNq
  echo $cubNq

  let Np=$Nq*$Nq*$Nq
  
  let maxE=6000000/$Np

  let skipE=$maxE/100
  
  echo $maxE
  
  for E in `seq 80 $skipE $maxE`
  do
    ./advectionMassMatrixMultiply $E
  done
done
      
