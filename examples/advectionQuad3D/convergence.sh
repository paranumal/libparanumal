#!/bin/bash

cd /scratch/stimmel/convergence
for meshsize in `seq 4 4 32`;
do
    for N in `seq 2 7`;
    do
	echo N=$N meshsize=$meshsize;
	~/holmes/examples/advectionQuad3D/advectionMainQuad3D $meshsize $N | grep norm ;
    done;
done
