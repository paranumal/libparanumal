/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "massHex3D.h"


void timeAxOperator(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options){

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	mesh_t *mesh = solver->mesh;
	
	// sync processes
	mesh->device.finish();
	MPI_Barrier(MPI_COMM_WORLD);
	
	double tic = MPI_Wtime();
	double AxTime;
	
	int iterations = 10;
	
	occa::streamTag start = mesh->device.tagStream();
	
	
#if 1
	
	void massOperator3D(solver_t *solver, dfloat lambda,
	                    occa::memory &o_q, occa::memory &o_Aq, const char *options);
	                    
	// assume 1 mpi process
	for(int it=0;it<iterations;++it){
	
		massOperator3D(solver, lambda, o_r, o_x, options);
	}
#else
	// assume 1 mpi process
	for(int it=0;it<iterations;++it)
		solver->partialAxKernel(solver->NlocalGatherElements,
		                        solver->o_localGatherElementList,
		                        solver->o_gjGeo,
		                        solver->o_gjD,
		                        solver->o_gjI,
		                        lambda, o_r,
		                        solver->o_grad,
		                        o_x);
#endif
		                        
	occa::streamTag end = mesh->device.tagStream();
	
	mesh->device.finish();
	double toc = MPI_Wtime();
	
	double localElapsed = toc-tic;
	//  localElapsed = mesh->device.timeBetween(start, end);
	
	// time cudamemcpy for same amount of data movement
	int gjNp = mesh->gjNq*mesh->gjNq*mesh->gjNq;
	int Nbytes =((sizeof(dfloat)*(2*mesh->Np*2+ 7*gjNp)/2)); // use 1/2 because of load+store
	occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
	occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);
	
	mesh->device.finish();
	tic = MPI_Wtime();
	
	occa::streamTag startCopy = mesh->device.tagStream();
	for(int it=0;it<iterations;++it){
		o_bah.copyTo(o_foo);
	}
	occa::streamTag endCopy = mesh->device.tagStream();
	
	mesh->device.finish();
	toc = MPI_Wtime();
	double copyElapsed = (toc-tic);
	copyElapsed = mesh->device.timeBetween(startCopy, endCopy);
	double copyBandwidth = mesh->Nelements*((Nbytes*iterations*2.)/(1024.*1024.*1024.*copyElapsed));
	
	
	int   localDofs = mesh->Np*mesh->Nelements;
	int localElements = mesh->Nelements;
	double globalElapsed;
	int   globalDofs;
	int   globalElements;
	int    root = 0;
	
	MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD );
	MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_INT,   MPI_SUM, root, MPI_COMM_WORLD );
	MPI_Reduce(&localElements,&globalElements,1, MPI_INT,   MPI_SUM, root, MPI_COMM_WORLD );
	
	int gjNq = mesh->gjNq;
	gjNp = mesh->gjNq*mesh->gjNq*mesh->gjNq;
	int Nq = mesh->Nq;
	
	double flops;
	double bw;
	if(!strstr(options, "COLLOCATION")){
		flops =
		  gjNq*Nq*Nq*Nq*2 +
		  gjNq*gjNq*Nq*Nq*2 +
		  gjNq*gjNq*gjNq*Nq*2 +
		  gjNq*gjNq*gjNq +
		  gjNq*gjNq*gjNq*Nq*2 +
		  gjNq*gjNq*Nq*Nq*2 +
		  gjNq*Nq*Nq*Nq*2 ;
		bw = sizeof(dfloat)*(2*Nq*Nq*Nq + gjNp);
	}else{
		printf("collocation\n");
		flops =
		  Nq*Nq*Nq*Nq*12 +
		  Nq*Nq*Nq*15;
		bw = sizeof(dfloat)*(3*Nq*Nq*Nq);
	}
	
	double gflops = globalElements*flops*iterations/(1024*1024*1024.*globalElapsed);
	bw *= (globalElements*iterations)/(1024.*1024*1024*globalElapsed);
	if(rank==root){
		printf("%02d %02d %02d %17.15lg %d %17.15E %17.15E %17.15E %17.15E %17.15E\t"
		       "[ RANKS N DOFS ELAPSEDTIME ITERATIONS (DOFS/RANKS) (DOFS/TIME/ITERATIONS/RANKS) (Ax GFLOPS) (copy GB/s) (achieved GB/s)]\n",
		       size, mesh->N, globalDofs, globalElapsed, iterations, globalDofs/(double)size,
		       (globalDofs*iterations)/(globalElapsed*size), gflops, copyBandwidth, bw );
	
printf("NUMBER OF INTEREST: %d %17.15E \n", mesh->Nq-1, gflops);
printf("COPYBW %d % 17.15E \n", mesh->Nq-1,copyBandwidth );
}
	
	
}

void timeSolver(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options){

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	mesh_t *mesh = solver->mesh;
	
	// sync processes
	mesh->device.finish();
	MPI_Barrier(MPI_COMM_WORLD);
	
	double tic = MPI_Wtime();
	int maxIterations = 3000;
	double AxTime;
	
	int iterations = massSolveHex3D(solver, lambda, o_r, o_x, maxIterations, options);
	
	mesh->device.finish();
	double toc = MPI_Wtime();
	
	double localElapsed = toc-tic;
	int   localDofs = mesh->Np*mesh->Nelements;
	double globalElapsed;
	int   globalDofs;
	int    root = 0;
	
	MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD );
	MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_INT,   MPI_SUM, root, MPI_COMM_WORLD );
	
	if(rank==root){
		printf("%02d %02d %02d %17.15lg %d %17.15E %17.15E \t [ RANKS N DOFS ELAPSEDTIME ITERATIONS (DOFS/RANKS) (DOFS/TIME/ITERATIONS/RANKS) \n",
		       size, mesh->N, globalDofs, globalElapsed, iterations, globalDofs/(double)size, (globalDofs*(double)iterations)/(globalElapsed*size));
	}
	
}




int main(int argc, char **argv){

	// start up MPI
	MPI_Init(&argc, &argv);
	
	if(argc!=3){
		// to run cavity test case with degree N elements
		printf("usage: ./main meshes/cavityH005.msh N\n");
		exit(-1);
	}
	
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	
	// int specify polynomial degree
	int N = atoi(argv[2]);
	
	// solver can be CG or PCG
	// preconditioner can be JACOBI, OAS, NONE
	// method can be CONTINUOUS or IPDG
	char *options =
	  strdup("solver=CG method=CONTINUOUS preconditioner=NONE");
	//strdup("solver=CG method=IPDG preconditioner=NONE");
	
	// set up mesh stuff
	mesh3D *mesh = meshSetupHex3D(argv[1], N);
	
	ogs_t *ogs;
	precon_t *precon;
	
	// parameter for mass problem (-laplacian + lambda)*q = f
	dfloat lambda = 1;
	
	// set up
	occa::kernelInfo kernelInfo;
	massSetupHex3D(mesh, kernelInfo);
	
	solver_t *solver = massSolveSetupHex3D(mesh, lambda, kernelInfo, options);
	int Nall = (mesh->Nq+1)*(mesh->Nq+1)*(mesh->Nq+1)*(mesh->Nelements+mesh->totalHaloPairs);
	//int Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
	printf("Nall = %d, mesh->Np = %d mesh->Nelements = %d mesh->totalHaloPairs = %d \n", Nall, mesh->Np, mesh->Nelements, mesh->totalHaloPairs);
	dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
	dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));
	
	// load rhs into r
	for(int e=0;e<mesh->Nelements;++e){
		for(int n=0;n<mesh->Np;++n){
		
			int ggid = e*mesh->Np*mesh->Nggeo + n;
			dfloat wJ = mesh->ggeo[ggid+mesh->Np*GWJID];
			
			int   id = e*mesh->Np+n;
			dfloat xn = mesh->x[id];
			dfloat yn = mesh->y[id];
			dfloat zn = mesh->z[id];
			
			dfloat f = cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);
			
			r[id] = wJ*f;
			
			x[id] = 0; // initial guess
		}
	}
	
	occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
	occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);
	
	timeAxOperator(solver, lambda, o_r, o_x, options);
	
	//  timeSolver(solver, lambda, o_r, o_x, options);
	
	// copy solution from DEVICE to HOST
	o_x.copyTo(mesh->q);
	
	dfloat maxError = 0;
	for(int e=0;e<mesh->Nelements;++e){
		for(int n=0;n<mesh->Np;++n){
			int   id = e*mesh->Np+n;
			dfloat xn = mesh->x[id];
			dfloat yn = mesh->y[id];
			dfloat zn = mesh->z[id];
			//printf("xd = %lf yn = %lf zn = %lf \n",xn,yn,zn);
			dfloat exact = cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);
			dfloat error = fabs(exact-mesh->q[id]);
#if 0
			if (error > 5){
				printf("element %d id = %d exact %lf comp. %lf error %lf \n",e, id, exact, mesh->q[id], exact-mesh->q[id]);
			}
#endif
			maxError = mymax(maxError, error);
			
			//mesh->q[id] -= exact;
		}
	}
	
	dfloat globalMaxError = 0.0f;
	MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
	if(rank==0)
		printf("globalMaxError = %17.15g\n", globalMaxError);
		
		
	// close down MPI
	MPI_Finalize();
	
	exit(0);
	return 0;
}
