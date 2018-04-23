#!/bin/bash

for N in `seq 70 128`;
do  
    ./cubed_grid $N > ../../../meshes/cubed_grid_${N}.msh
    echo "$N complete"
done
