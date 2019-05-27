#!/bin/bash

rm massMatrixMultiplyVT
make massMatrixMultiplyVT

# to run a sweep
for Nq in `seq 2 1 11`
do

    #  let cubNq=$(($Nq + 2))
  let cubNq=$((1+$Nq))
  echo $cubNq

  let Np=$Nq*$Nq*$Nq
  
  let maxE=2000000/$Np

  let tmpE=$maxE/400
  let tmpE2=(19+$tmpE)/20
  let skipE=$tmpE2*20
  
  echo $maxE
  echo $skipE

  ./massMatrixMultiplyVT $Nq $cubNq 1  1
  for mode in `seq 1 4`
  do
      for E in `seq 80 $skipE $maxE`
      do
	  ./massMatrixMultiplyVT $Nq $cubNq $E $mode
      done
  done
done
      
