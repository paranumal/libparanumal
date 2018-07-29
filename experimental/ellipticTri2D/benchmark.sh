#!/bin/bash

#make clean; make -j

declare -a testName=(
"CG                           "
"CG+Jacobi                    " 
"AMG,vcycle                   " 
"AMG,vcycle,Cheb              "
"AMG,kcycle                   " 
"AMG,kcycle,Cheb              "  
"OAS,Jacobi                   "
"OAS,EXACTFULLPATCH           "
"OAS,APPROXFULLPATCH          "
"OAS,EXACTFACE                "
"OAS,APPROXFACE               "
"OAS,EXACTBLOCK               "
"OAS,APPROXBLOCK              "
"HybAMG,DAMPEDJACOBI          "
"HybAMG,DAMPEDJACOBI,CHEB     "
"HybAMG,EXACTFULLPATCH        "
"HybAMG,EXACTFULLPATCH,CHEB   "
"HybAMG,APPROXFULLPATCH       "
"HybAMG,APPROXFULLPATCH,CHEB  "
"HybAMG,EXACTFACEPATCH        "
"HybAMG,EXACTFACEPATCH,CHEB   "
"HybAMG,APPROXFACEPATCH       "
"HybAMG,APPROXFACEPATCH,CHEB  "
"HybAMG,EXACTBLOCKJACOBI      "
"HybAMG,EXACTBLOCKJACOBI,CHEB "
"HybAMG,APPROXBLOCKJACOBI     "
"HybAMG,APPROXBLOCKJACOBI,CHEB"
)

declare -a options=(
"solver=PCG,FLEXIBLE method=IPDG preconditioner=NONE"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=JACOBI" 

"solver=PCG,FLEXIBLE method=IPDG preconditioner=FULLALMOND" 
"solver=PCG,FLEXIBLE method=IPDG preconditioner=FULLALMOND"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=FULLALMOND" 
"solver=PCG,FLEXIBLE method=IPDG preconditioner=FULLALMOND"  

"solver=PCG,FLEXIBLE method=IPDG preconditioner=OAS smoother=DAMPEDJACOBI"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=OAS smoother=EXACTFULLPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=OAS smoother=APPROXFULLPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=OAS smoother=EXACTFACEPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=OAS smoother=APPROXFACEPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=OAS smoother=EXACTBLOCKJACOBI"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=OAS smoother=APPROXBLOCKJACOBI"

"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=DAMPEDJACOBI"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=DAMPEDJACOBI,CHEBYSHEV"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=EXACTFULLPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=EXACTFULLPATCH,CHEBYSHEV"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=APPROXFULLPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=APPROXFULLPATCH,CHEBYSHEV"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=EXACTFACEPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=EXACTFACEPATCH,CHEBYSHEV"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=APPROXFACEPATCH"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=APPROXFACEPATCH,CHEBYSHEV"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=EXACTBLOCKJACOBI"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=EXACTBLOCKJACOBI,CHEBYSHEV"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=APPROXBLOCKJACOBI"
"solver=PCG,FLEXIBLE method=IPDG preconditioner=MULTIGRID,HALFDOFS smoother=APPROXBLOCKJACOBI,CHEBYSHEV"
)

declare -a parAlmondOptions=(
"solver= smoother= partition=" 
"solver= smoother= partition=" 

"solver=VCYCLE smoother=DAMPEDJACOBI partition=STRONGNODES" 
"solver=VCYCLE smoother=CHEBYSHEV    partition=STRONGNODES" 
"solver=KCYCLE smoother=DAMPEDJACOBI partition=STRONGNODES" 
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"

"solver=EXACT  smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=EXACT  smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=EXACT  smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=EXACT  smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=EXACT  smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=EXACT  smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=EXACT  smoother=CHEBYSHEV    partition=STRONGNODES"

"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
"solver=KCYCLE smoother=CHEBYSHEV    partition=STRONGNODES"
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
