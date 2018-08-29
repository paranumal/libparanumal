#!/bin/bash

for N in `seq 3 8`;
do
    for alpha in `seq 0 10`;
    do
	echo N=$N alpha=$alpha;
	./advectionMainQuad3D 20 $N $alpha | grep norm ;
    done;
done
