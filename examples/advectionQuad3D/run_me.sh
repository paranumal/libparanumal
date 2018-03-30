#!/bin/bash

#Currently runs lserk only, filter enabled, alpha = 1./N, with normal upwind coefficients
#Should give ~1e-9 when stable, ~1e-2 filter only, ~1e9 with alpha=0
#Intermediate values are plotted every 1000 steps

./advectionMainQuad3D /home/stimmel/holmes/meshes/cubed_grid_28.msh 5
