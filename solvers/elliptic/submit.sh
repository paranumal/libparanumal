# LSF directives                                                                                                                                                               

#BSUB -P CSC262                                                                                                                                                                                                                                                                                                                                           
#BSUB -J paranumal                                                                                                                                                                                                                                                                                                                                            
#BSUB -o bp5.o%J                                                                                                                                                                                                                                                                                                                                              
#BSUB -W 1:00

mesh_file=$1
N=$2
maxiter=$3
nodes=$4


cd  /ccs/home/karakus/libparanumal/solvers/elliptic

# 1 First run this script to generate  OCCA kernels for only once
./BP  ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}

# 2.After OCCA kernels are built, following script can be  uses
#jsrun -n${nodes} -r1 -a4 -c4 -g4 ./BP  ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}


#jsrun -n${nodes} -r1 -a2 -c2 -g2 ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter                                                                                       
                                                                                                                                                                               
#jsrun -n${nodes} -r1 -a1 -c1 -g1 ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter                                                                                     
                                                                                                                                                                              
#jsrun -n${nodes} -r1 -a1 -c1 -g1 nvprof --log-file bp5.${nodes}.%p.log --profile-from-start-off ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter                        
                                                                                                                                                                              
#jsrun -n${nodes} -r1 -a2 -c2 -g2 --smpiargs "-gpu" ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter                                                                                                                                                                                                                                                    
#jsrun -n${nodes} -r1 -a4 -c4 -g4 --smpiargs "-gpu" ./BPMain ./setups/setupHex3D.rc $mesh_file $N $maxiter

#jsrun -n${nodes} -r1 -a4 -c4 -g4 nvprof --annotate-mpi openmpi -o BP5_N1.%q{OMPI_COMM_WORLD_RANK}.prof  ./BP ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}

#jsrun -n${nodes} -r1 -a4 -c4 -g4 ./BP  ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}

#jsrun -n${nodes} -r1 -a4 -c4 -g4 ./ellipticMain  ./setups/setupHex3D.rc 


#jsrun -n${nodes} -r1 -a4 -c4 -g4 nvprof --annotate-mpi openmpi -o BP5_N1.%q{OMPI_COMM_WORLD_RANK}.prof  ./BP ./setups/setupHex3D.rc ${mesh_file} ${N} ${maxiter}
#./ellipticMain  ./setups/setupHex3D.rc 
