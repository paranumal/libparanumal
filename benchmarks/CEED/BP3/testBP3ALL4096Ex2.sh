#! /bin/bash
#
#PBS -l walltime=00:25:00
#PBS -l nodes=1:ppn=28:gpus=1
#PBS -W group_list=newriver
#PBS -A p100_test
#PBS -q p100_normal_q

#PBS -j oe

cd $PBS_O_WORKDIR

module purge
module load gcc/5.2.0 cuda openmpi 
sed '372s/.*/"ellipticAxHex3D_Ref3D0",/' ellipticSolveSetupHex3D.c > foo.c
mv foo.c ellipticSolveSetupHex3D.c
make -j

for n in {1,2,3,4,5,6,7,8} 
do    
./ellipticMainHex3D ../../../meshes/cubeHexH0125.msh $n  
done

sed '372s/.*/"ellipticAxHex3D_Ref3D1",/' ellipticSolveSetupHex3D.c > foo.c
mv foo.c ellipticSolveSetupHex3D.c
make -j

for n in {1,2,3,4,5,6,7,8} 
do
./ellipticMainHex3D ../../../meshes/cubeHexH0125.msh $n
done

sed '372s/.*/"ellipticAxHex3D_Ref3D2",/' ellipticSolveSetupHex3D.c > foo.c
mv foo.c ellipticSolveSetupHex3D.c
make -j

for n in {1,2,3,4,5,6,7,8} 
do
./ellipticMainHex3D ../../../meshes/cubeHexH0125.msh $n
done

sed '372s/.*/"ellipticAxHex3D_Ref3D3",/' ellipticSolveSetupHex3D.c > foo.c
mv foo.c ellipticSolveSetupHex3D.c
make -j

for n in {1,2,3,4,5,6,7,8} 
do
./ellipticMainHex3D ../../../meshes/cubeHexH0125.msh $n
done

sed '372s/.*/"ellipticAxHex3D_Ref3D4",/' ellipticSolveSetupHex3D.c > foo.c
mv foo.c ellipticSolveSetupHex3D.c
make -j

for n in {1,2,3,4,5,6,7,8} 
do
./ellipticMainHex3D ../../../meshes/cubeHexH0125.msh $n
done

sed '372s/.*/"ellipticAxHex3D_Ref3D5",/' ellipticSolveSetupHex3D.c > foo.c
mv foo.c ellipticSolveSetupHex3D.c
make -j

for n in {1,2,3,4,5,6,7,8} 
do
./ellipticMainHex3D ../../../meshes/cubeHexH0125.msh $n
done

sed '372s/.*/"ellipticAxHex3D_Ref3D6",/' ellipticSolveSetupHex3D.c > foo.c
mv foo.c ellipticSolveSetupHex3D.c
make -j

for n in {1,2,3,4,5,6,7,8} 
do
./ellipticMainHex3D ../../../meshes/cubeHexH0125.msh $n
done

sed '372s/.*/"ellipticAxHex3D_Ref3D7",/' ellipticSolveSetupHex3D.c > foo.c
mv foo.c ellipticSolveSetupHex3D.c
make -j

for n in {1,2,3,4,5,6,7,8} 
do
./ellipticMainHex3D ../../../meshes/cubeHexH0125.msh $n
done

sed '372s/.*/"ellipticAxHex3D_Ref3D8",/' ellipticSolveSetupHex3D.c > foo.c
mv foo.c ellipticSolveSetupHex3D.c
make -j

for n in {1,2,3,4,5,6,7,8} 
do
./ellipticMainHex3D ../../../meshes/cubeHexH0125.msh $n
done



