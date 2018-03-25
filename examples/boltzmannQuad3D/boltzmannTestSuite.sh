#!/bin/bash

cd /scratch/stimmel/mrab
export WORK_DIR="/home/stimmel/holmes/examples/boltzmannQuad3D"
export MESH_DIR="/home/stimmel/holmes/meshes"
#$WORK_DIR/boltzmannMainQuad3D $MESH_DIR/cubed_grid_layers_medium.msh 2 1 0 0 1 1
#$WORK_DIR/boltzmannMainQuad3D $MESH_DIR/cubed_grid_layers_medium.msh 2 0 1 0 1 1
$WORK_DIR/boltzmannMainQuad3D $MESH_DIR/cubed_grid_layers_medium.msh 2 1 1 0 1 1
#$WORK_DIR/boltzmannMainQuad3D $MESH_DIR/cubed_grid_layers_medium.msh 2 1 1 0 1 2
#$WORK_DIR/boltzmannMainQuad3D $MESH_DIR/cubed_grid_layers_medium.msh 2 1 1 1 1 2
#$WORK_DIR/boltzmannMainQuad3D $MESH_DIR/cubed_grid_layers_medium.msh 2 1 1 1 2 2
