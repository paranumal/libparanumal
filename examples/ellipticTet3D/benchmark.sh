#!/bin/bash

#make clean; make -j

declare -a testName=(
"CG+Jacobi                    "

"HybAMG,DAMPEDJACOBI          "
"HybAMG,DAMPEDJACOBI,CHEB     "
"HybAMG,APPROXFULLPATCH       "
"HybAMG,APPROXFACEPATCH       "
"HybAMG,APPROXBLOCKJACOBI     "
)

declare -a options=(
"solver=PCG,FLEXIBLE method=IPDG preconditioner=JACOBI"

"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=DAMPEDJACOBI"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=DAMPEDJACOBI,CHEBYSHEV"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=APPROXFULLPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=APPROXFACEPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=APPROXBLOCKJACOBI"
)

declare -a parAlmondOptions=(
"solver= smoother= partition="

"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
)

numTests=${#options[@]}

declare -a meshes=(
"../../meshes/cavity3DH05.msh"
"../../meshes/cavity3DH01.msh"
"../../meshes/cavity3DH005.msh"
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
