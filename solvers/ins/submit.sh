# LSF directives                                                                                                                                                               

#BSUB -P CSC262                                                                                                                                                                                                                                                                                                                                           
#BSUB -J paranumal                                                                                                                                                                                                                                                                                                                                               
#BSUB -o ins.o%J                                                                                                                                                                                                                                                                                                                                              
#BSUB -W 1:00

mesh_file=$1
N=$2
maxiter=$3
nodes=$4

# End LSF directives and begin shell commands                                                                                                                                                                                                                                                                                                                 
cd  /ccs/home/karakus/libparanumal/solvers/ins
#jsrun -n${nodes} -r1 -a2 -c2 -g2 ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter                                                                                       
                                                                                                                                                                               
#jsrun -n${nodes} -r1 -a1 -c1 -g1 ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter                                                                                     
                                                                                                                                                                              
#jsrun -n${nodes} -r1 -a1 -c1 -g1 nvprof --log-file bp5.${nodes}.%p.log --profile-from-start-off ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter                        
                                                                                                                                                                              
#jsrun -n${nodes} -r1 -a2 -c2 -g2 --smpiargs "-gpu" ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter                                                                                                                                                                                                                                                    
#jsrun -n${nodes} -r1 -a4 -c4 -g4 --smpiargs "-gpu" ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter

#jsrun -n${nodes} -r1 -a4 -c4 -g4 nvprof --annotate-mpi openmpi -o BP5_N1.%q{OMPI_COMM_WORLD_RANK}.prof  ./BP ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}

#jsrun -n${nodes} -r1 -a4 -c4 -g4 nvprof --annotate-mpi openmpi -o BP5_N1.%q{OMPI_COMM_WORLD_RANK}.prof  ./BP ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}
jsrun -n${nodes} -r1 -a4 -c4 -g4  ./insTest ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}

#./insTest setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}
#./insMain ./setups/setupHex3D.rc
