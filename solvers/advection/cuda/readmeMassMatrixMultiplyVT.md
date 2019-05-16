#!/bin/bash

rm massMatrixMultiplyVT
make massMatrixMultiplyVT

# to run a sweep
for Nq in `seq 8 2 12`
do

  let cubNq=$(($Nq + 2))
  echo $cubNq

  let Np=$Nq*$Nq*$Nq
  
  let maxE=2000000/$Np

  let tmpE=$maxE/400
  let tmpE2=(19+$tmpE)/20
  let skipE=$tmpE2*20
  
  echo $maxE
  echo $skipE

  ./massMatrixMultiplyVT $Nq $cubNq 1 
  
  for E in `seq 80 $skipE $maxE`
  do
    ./massMatrixMultiplyVT $Nq $cubNq $E
  done
done
      
