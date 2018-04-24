#!/bin/bash

cd /scratch/stimmel/convergence
for N in `seq 2 9`;
do
    for alpha in `seq 0 1 10`;
    do
	echo N=$N alpha=$alpha;
	~/holmes/examples/advectionQuad3D/advectionMainQuad3D ~/holmes/meshes/cubed_grid_28.msh $N $alpha | grep norm ;
    done;
done
