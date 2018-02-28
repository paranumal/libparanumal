#!/bin/bash

#make clean; make -j

declare -a testName=( 
"AMG,kcycle,Cheb              "  
"SEMFEM,AMG,kcycle,Cheb       "  
)

declare -a options=(
"solver=PCG,FLEXIBLE method=CONTINUOUS preconditioner=FULLALMOND" 
"solver=PCG,FLEXIBLE method=CONTINUOUS preconditioner=SEMFEM" 
)

declare -a parAlmondOptions=(
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES" 
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
)

numTests=${#options[@]}

declare -a meshes=(
"../../meshes/cavityH05.msh"
"../../meshes/cavityH01.msh"
"../../meshes/cavityH005.msh"
"../../meshes/cavityH0025.msh"
"../../meshes/cavityH00125.msh"
)

numMeshes=${#meshes[@]}

for (( i=1; i<${numTests}+1; i++ ));
do
  for (( m=1; m<${numMeshes}+1; m++ ));
  do
    echo ${meshes[$m-1]}
    for (( p=1; p<9; p++ ));
    do
      echo "Test " $i " : " ${testName[$i-1]}
      ./ellipticMainTri2D ${meshes[$m-1]} $p homogeneous2D.h "${options[$i-1]}" "${parAlmondOptions[$i-1]}" |& grep RAN 
    done
  done
done
