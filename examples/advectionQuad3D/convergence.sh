#!/bin/bash

for N in `seq 3 8`;
do
    for h in `seq 4 4 28`;
    do
	echo N=$N h=$h;
	./advectionMainQuad3D $h $N 1 | grep norm ;
    done;
done
