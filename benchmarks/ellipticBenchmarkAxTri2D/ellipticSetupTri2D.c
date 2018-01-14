#include "ellipticBenchmarkTri2D.h"

solver_t *ellipticSetupTri2D(mesh_t *mesh, dfloat tau, dfloat lambda, iint*BCType,
		occa::kernelInfo &kernelInfo, const char *options, const char *parAlmondOptions,
		int Nblocks, int Nnodes){

	mesh->Nfields = 1;

	// OCCA build stuff
	char deviceConfig[BUFSIZ];
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// use rank to choose DEVICE
	sprintf(deviceConfig, "mode = CUDA, deviceID = %d", rank%2);
	//sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
	//sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 1);
	//sprintf(deviceConfig, "mode = Serial");

	meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

	iint Ntotal = mesh->Np*mesh->Nelements;
	iint Nblock = (Ntotal+blockSize-1)/blockSize;
	iint Nhalo = mesh->Np*mesh->totalHaloPairs;
	iint Nall   = Ntotal + Nhalo;

	solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));

	solver->tau = tau;

	solver->mesh = mesh;

	solver->p   = (dfloat*) calloc(Nall,   sizeof(dfloat));
	solver->z   = (dfloat*) calloc(Nall,   sizeof(dfloat));
	solver->Ax  = (dfloat*) calloc(Nall,   sizeof(dfloat));
	solver->Ap  = (dfloat*) calloc(Nall,   sizeof(dfloat));
	solver->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));

	solver->grad = (dfloat*) calloc(Nall*4, sizeof(dfloat));

	solver->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
	solver->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
	solver->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);

	solver->o_res = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
	solver->o_Sres = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
	solver->o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
	solver->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), solver->Ap);
	solver->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), solver->tmp);

	solver->o_grad  = mesh->device.malloc(Nall*4*sizeof(dfloat), solver->grad);

	//setup async halo stream
	solver->defaultStream = mesh->device.getStream();
	solver->dataStream = mesh->device.createStream();
	mesh->device.setStream(solver->defaultStream);

	iint Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
	if(Nbytes>0){
		occa::memory o_sendBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
		occa::memory o_recvBuffer = mesh->device.mappedAlloc(Nbytes, NULL);

		solver->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
		solver->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();

		occa::memory o_gradSendBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);
		occa::memory o_gradRecvBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);

		solver->gradSendBuffer = (dfloat*) o_gradSendBuffer.getMappedPointer();
		solver->gradRecvBuffer = (dfloat*) o_gradRecvBuffer.getMappedPointer();
	}else{
		solver->sendBuffer = NULL;
		solver->recvBuffer = NULL;
	}
	mesh->device.setStream(solver->defaultStream);

	solver->type = strdup(dfloatString);

	solver->Nblock = Nblock;

	//fill geometric factors in halo
	if(mesh->totalHaloPairs){
		iint Nlocal = mesh->Nelements*mesh->Np;
		iint Nhalo  = mesh->totalHaloPairs*mesh->Np;
		dfloat *vgeoSendBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Nvgeo, sizeof(dfloat));

		// import geometric factors from halo elements
		mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nvgeo*sizeof(dfloat));

		meshHaloExchange(mesh,
				mesh->Nvgeo*sizeof(dfloat),
				mesh->vgeo,
				vgeoSendBuffer,
				mesh->vgeo + mesh->Nelements*mesh->Nvgeo);

		mesh->o_vgeo =
			mesh->device.malloc((mesh->Nelements + mesh->totalHaloPairs)*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);
	}

	//check all the bounaries for a Dirichlet
	bool allNeumann = (lambda==0) ? true :false;
	solver->allNeumannPenalty = 1;
	iint totalElements = 0;
	MPI_Allreduce(&(mesh->Nelements), &totalElements, 1, MPI_IINT, MPI_SUM, MPI_COMM_WORLD);
	solver->allNeumannScale = 1.0/sqrt(mesh->Np*totalElements);

	solver->EToB = (int *) calloc(mesh->Nelements*mesh->Nfaces,sizeof(int));
	for (iint e=0;e<mesh->Nelements;e++) {
		for (int f=0;f<mesh->Nfaces;f++) {
			int bc = mesh->EToB[e*mesh->Nfaces+f];
			if (bc>0) {
				int BC = BCType[bc]; //get the type of boundary
				solver->EToB[e*mesh->Nfaces+f] = BC; //record it
				if (BC!=2) allNeumann = false; //check if its a Dirchlet
			}
		}
	}
	MPI_Allreduce(&allNeumann, &(solver->allNeumann), 1, MPI::BOOL, MPI_LAND, MPI_COMM_WORLD);
	printf("allNeumann = %d \n", solver->allNeumann);

	solver->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int), solver->EToB);

	//add standard boundary functions
	char *boundaryHeaderFileName;
	boundaryHeaderFileName = strdup(DHOLMES "/examples/ellipticTri2D/ellipticBoundary2D.h");
	kernelInfo.addInclude(boundaryHeaderFileName);

	kernelInfo.addParserFlag("automate-add-barriers", "disabled");

	if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
		kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
	}

	kernelInfo.addDefine("p_blockSize", blockSize);

	// add custom defines
	kernelInfo.addDefine("p_NpP", (mesh->Np+mesh->Nfp*mesh->Nfaces));
	kernelInfo.addDefine("p_Nverts", mesh->Nverts);
	//sizes for the coarsen and prolongation kernels. degree N to degree 1
	kernelInfo.addDefine("p_NpFine", mesh->Np);
	kernelInfo.addDefine("p_NpCoarse", mesh->Nverts);

	int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
	kernelInfo.addDefine("p_Nmax", Nmax);

	int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
	kernelInfo.addDefine("p_maxNodes", maxNodes);

	kernelInfo.addDefine("p_NblockV", Nblocks);
	kernelInfo.addDefine("p_NnodesV", Nnodes);
	kernelInfo.addDefine("p_maxNnzPerRow", mesh->maxNnzPerRow);

	kernelInfo.addDefine("p_NblockS", Nblocks);
	kernelInfo.addDefine("p_NblockP", Nblocks);
	kernelInfo.addDefine("p_NblockG", Nblocks);
	iint *globalGatherElementList    = (iint*) calloc(mesh->Nelements, sizeof(iint));
	iint *localGatherElementList = (iint*) calloc(mesh->Nelements, sizeof(iint));
	int globalCount = 0;
	int localCount =0;

	for(iint e=0;e<mesh->Nelements;++e){
		globalGatherElementList[globalCount++] = e;
		localGatherElementList[localCount++] = e;
	}

	solver->NglobalGatherElements = mesh->Nelements;
	solver->NlocalGatherElements = mesh->Nelements;

	solver->o_globalGatherElementList =
		mesh->device.malloc(globalCount*sizeof(iint), globalGatherElementList);

	solver->o_localGatherElementList =
		mesh->device.malloc(localCount*sizeof(iint), localGatherElementList);
	mesh->invSparseV = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
	for (int n=0;n<mesh->Np*mesh->Np;n++)
		mesh->invSparseV[n] = mesh->sparseV[n];

	int paddedRowSize = 4*((mesh->SparseNnzPerRow+3)/4); //make the nnz per row a multiple of 4

	char* IndTchar = (char*) calloc(paddedRowSize*mesh->Np,sizeof(char));
	for (int m=0;m<paddedRowSize/4;m++) {
		for (int n=0;n<mesh->Np;n++) {
			for (int k=0;k<4;k++) {
				if (k+4*m < mesh->SparseNnzPerRow) {
					IndTchar[k+4*n+m*4*mesh->Np] = mesh->sparseStackedNZ[n+(k+4*m)*mesh->Np];
				} else {
					IndTchar[k+4*n+m*4*mesh->Np] = 0;
				}
			}
		}
	}



	kernelInfo.addDefine("p_SparseNnzPerRow", paddedRowSize);
	printf("sparse NNZ %d \n", paddedRowSize);  


	mesh->o_sparseSrrT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrrT);
	mesh->o_sparseSrsT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSrsT);
	mesh->o_sparseSssT = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(dfloat), mesh->sparseSssT);

	mesh->SparseNnzPerRow = paddedRowSize;
	mesh->o_sparseStackedNZ = mesh->device.malloc(mesh->Np*mesh->SparseNnzPerRow*sizeof(char), IndTchar);
	free(IndTchar);



	return solver;
}


