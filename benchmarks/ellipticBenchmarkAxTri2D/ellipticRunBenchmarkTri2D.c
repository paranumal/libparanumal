#include "ellipticBenchmarkTri2D.h"

void ellipticRunBenchmark2D(solver_t *solver, char *options, occa::kernelInfo kernelInfo, char *kernelFileName, int Nblocks, int Nnodes){

	mesh2D *mesh = solver->mesh;

	int Ntrials = 10;

	size_t L1CacheSize = 24576; //L1 cache size of test device in bytes (24KB for P100)

	int NKernels;
	char kernelName[BUFSIZ];

	NKernels = 1;
	sprintf(kernelName, "ellipticPartialAxTri2D");

	//  kernelInfo.addCompilerFlag("-G");

	dfloat time = 0.;

	char testkernelName[BUFSIZ];
	occa::kernel testKernel;
	for(iint i=0; i<NKernels; i++) {

		sprintf(testkernelName, "%s_v%d", kernelName,  i);
		printf("%s Kernel #%02d\n", kernelFileName, i);
		printf("Np = %d sizeof(dfloat) = %d Nelements = %d \n", mesh->Np, sizeof(dfloat), mesh->Nelements);

		testKernel = mesh->device.buildKernelFromSource(kernelFileName,testkernelName,kernelInfo);

		dfloat lambda = 0;
		// sync processes

		iint Nbytes = sizeof(dfloat)*mesh->Np*2; // load one field, save one filed (ignore geofacs)
		Nbytes /= 2;
		occa::memory o_foo = mesh->device.malloc(Nbytes*mesh->Nelements);
		occa::memory o_bah = mesh->device.malloc(Nbytes*mesh->Nelements);
		mesh->device.finish();

		occa::streamTag startCopy = mesh->device.tagStream();
		for(int it=0;it<Ntrials;++it){
			o_bah.copyTo(o_foo);
		}
		occa::streamTag endCopy = mesh->device.tagStream();
		mesh->device.finish();

		double  copyElapsed = mesh->device.timeBetween(startCopy, endCopy);
		printf("copied \n");

		// time Ax
		double timeAx = 0.0f;
		double kernelElapsed =0.0f;;
		occa::streamTag start[Ntrials], end[Ntrials];
		for(int it=0;it<Ntrials;++it){
			mesh->device.finish();
			start[it] = mesh->device.tagStream();

			testKernel(mesh->Nelements, 
					solver->o_localGatherElementList,
					mesh->o_ggeo, 
					mesh->o_sparseStackedNZ, 
					mesh->o_sparseSrrT, 
					mesh->o_sparseSrsT, 
					mesh->o_sparseSssT,		 
					mesh->o_MM, 
					lambda,
					solver->o_p,
					solver->o_grad);


			end[it] = mesh->device.tagStream();
			timeAx = mesh->device.timeBetween(start[it],end[it]);
			kernelElapsed +=timeAx;    
		}
		mesh->device.finish();
		kernelElapsed = kernelElapsed/(dfloat)Ntrials;

		//start
		//
		iint size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		iint   localDofs = mesh->Np*mesh->Nelements;
		iint   localElements = mesh->Nelements;
		iint   globalDofs;
		iint   globalElements;
		double globalCopyElapsed;

		MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_IINT,   MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce(&localElements,&globalElements,1, MPI_IINT,   MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce(&copyElapsed,&globalCopyElapsed,1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD );


		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank==0){
			printf("Ntrials = %d global copy el %f local copy el %f size(dfloat) %d \n",Ntrials, globalCopyElapsed, copyElapsed, sizeof(dfloat));

			// count actual number of non-zeros
			int nnzs = 0;
			for(iint n=0;n<mesh->Np*mesh->SparseNnzPerRow;++n)
				nnzs += (mesh->sparseStackedNZ[n]>0);
			printf("nnzs = %d\n", nnzs);

			// 6 flops per non-zero plus chain rule
			iint flops = nnzs*6 + mesh->Np*5;

			double roofline = ((mesh->Nelements*flops*(double)Ntrials))/(1e9*globalCopyElapsed);
			printf("Nelements = %d flops = %d Ntrials = %d copy elapsed scaled = %f kernelElapsed %16.16f\n",mesh->Nelements, flops, Ntrials,  1e9*globalCopyElapsed, kernelElapsed);

			double copyBandwidth   = mesh->Nelements*((Nbytes*Ntrials*2.)/(1e9*globalCopyElapsed));
			double kernelBandwidth = mesh->Nelements*((Nbytes*2.)/(1e9*kernelElapsed));
			double kernelGFLOPS = mesh->Nelements*flops/(1e9*kernelElapsed);

			printf("Copy Bw %16.17g achieved BW %16.17g\n", copyBandwidth, kernelBandwidth);
			printf("ROOFLINE %16.17g \n", roofline);
			printf("GFLOPS %16.17f \n", kernelGFLOPS);
			printf("time per kernel %f \n",kernelElapsed);
			printf("PARAMETERS %d %d %16.17f \n ", Nblocks, Nnodes, kernelGFLOPS );


		}
	}
}
