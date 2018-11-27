#mpiexec -np 1 ./ellipticMain ./setups/setupHex3D.rc
#mpiexec -np 2 ./ellipticMain ./setups/setupHex3D.rc
#mpiexec -np 4 ./ellipticMain ./setups/setupHex3D.rc
#mpiexec -np 8 ./ellipticMain ./setups/setupHex3D.rc

maxiter=100;
mesh_base=../../meshes/cubeHex;
mesh_extension=.msh;
for size in 1 2 4 8; do
    for mesh_id in 10 12; do
	for N in 2 3 4 5 6 7 8 9 11 15 ;do
	    meshfile=${mesh_base}${mesh_id}${mesh_extension}
            echo Reading ${meshfile} 
	    mpiexec -np ${size} ./BP ./setups/setupHex3D.rc ${meshfile} ${N} ${maxiter};
done;done;done;
