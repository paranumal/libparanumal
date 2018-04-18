#!/bin/bash

for N in `seq 4 70`;
do  
    ./cubed_grid $N > ../../../meshes/cubed_grid_${N}.msh
done
