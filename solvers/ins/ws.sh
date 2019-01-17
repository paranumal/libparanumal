#!/bin/bash                                                                                                                                                                     

function log2() {
  local input=$1
  local output=0

  while ((input>1)); do
    ((input/=2, output +=1))
  done

  echo $output
}

#for i in `seq 1 1 16`; do                                                                                                                                                      
#  a=$(log2 $i)                                                                                                                                                                 
#  echo "log of $i is $a"                                                                                                                                                       
#done                                                                                                                                                                           

maxiter=100;
mesh_base=../../meshes/box;
mesh_extension=.msh;
starting_elements=15


#for nodes in 1 2 4 8 16 32; do                                                                                                                                                 
#for nodes in 1 2 4 8 16 32; do                                                                                                                                                 
#for mesh_inc in 0; do                                                                
for nodes in 2; do
for N in 7; do
 # ((mesh_id=starting_elements+mesh_inc))
 # meshfile=${mesh_base}${mesh_id}${mesh_extension}
 # echo "Running ${meshfile} on ${nodes} node(s) at order ${N}"
  log_nodes=$(log2 $nodes)
  ((mesh_id=starting_elements+log_nodes))
  meshfile=${mesh_base}${mesh_id}${mesh_extension}
  echo "Running ${meshfile} on ${nodes} node(s) at order ${N}"
  bsub -nnodes ${nodes} -W 1:00 -P CSC262 -J paranumal -o insWeak.o%J \
        ./submit.sh ${meshfile} ${N} ${maxiter} ${nodes}
  sleep 3;
done;
done;
#done;
