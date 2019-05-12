#!/bin/bash

# to run a sweep
for Nq in `seq 2 2 12`
do

  rm advectionMassMatrixMultiply

  let cubNq=$(($Nq + 2))
  make advectionMassMatrixMultiply comp_Nq=$Nq comp_cubNq=$cubNq
  echo $cubNq

  let Np=$Nq*$Nq*$Nq
  
  let maxE=2000000/$Np

  let tmpE=$maxE/400
  let tmpE2=(19+$tmpE)/20
  let skipE=$tmpE2*20
  
  echo $maxE
  echo $skipE

  ./roofline "./advectionMassMatrixMultiply 1"
  
  for E in `seq 80 $skipE $maxE`
  do
    ./roofline "./advectionMassMatrixMultiply $E"
  done

  mv advectionMassMatrixMultiplyBlockedGlobalVgeoKernel.dat advectionMassMatrixMultiplyBlockedGlobalVgeoKernelNq${Nq}.dat
done
      
