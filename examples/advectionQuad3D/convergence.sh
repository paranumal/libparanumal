#!/bin/bash

for N in `seq 2 7`;
do
    for alpha in `seq 0 10`;
    do
	echo N=$N alpha=$alpha;
	./advectionMainQuad3D 20 $N $alpha | grep norm ;
    done;
done
