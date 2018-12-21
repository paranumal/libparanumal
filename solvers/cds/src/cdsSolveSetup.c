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

#include "cds.h"


void cdsSolveSetup(cds_t *cds, setupAide options){
  // sanity checking
  if (options.compareArgs("ELEMENT TYPE","HEXAHEDRA")){
      printf("ERROR: Only HEXAHEDRAL element is available currently\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
   }

  // Force cds to use the smae mesh with flow solver
   mesh_t *mesh      = cds->fSolver->mesh;
   
   ins_t  *fSolver   = cds->fSolver; 

   dlong Ntotal = mesh->Np*mesh->Nelements;
   dlong Nblock = mymax(1,(Ntotal+blockSize-1)/blockSize);
   dlong Nhalo  = mesh->Np*mesh->totalHaloPairs;
   dlong Nall   = Ntotal + Nhalo;

    cds->vOffset     = Ntotal; // keep it different for now.....b
    cds->sOffset     = Ntotal;
    // Properties have to be set from flow solver side
    cds->NSfields = 1; // hard coded right now.....
    // cds->NVfields = (mesh->dim == 2) ? 2 : 3; // get dimension from mesh
#if 0  // Use flow setting  
    cds->dt       = fSolver->dt;
    cds->g0       = fSolver->g0; 
    cds->NVfields = fSolver->NVfields;
    cds->rho      = fSolver->rho;
    cds->nu       = fSolver->nu; 
    //Reynolds number
    cds->Re       = fSolver->Re;          // Reynolds number

    // Get Solver parameters
    options.getArgs("CONDUCTIVITY", cds->k);
    options.getArgs("SPECIFIC HEAT", cds->cp);
    
    cds->alf       = cds->k/(cds->rho*cds->cp); // Thermal diffusivity
    cds->ialf      = 1.0/cds->alf;              // Inverse diff. 
    cds->idt       = 1.0/cds->dt;
    cds->lambda    = cds->g0 / (cds->dt * cds->alf);
    cds->Pr        = cds->nu/cds->alf;          // Prandtl number   

    // Assumes NStages are initiliazed in solver side.... 
    cds->S     = (dfloat*) calloc(cds->NSfields*cds->Nstages*Ntotal,sizeof(dfloat));
    //additional field storage
    cds->NS    = (dfloat*) calloc(cds->NSfields*(cds->Nstages+1)*Ntotal,sizeof(dfloat));
    cds->rkS   = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));

    if(cds->Nsubsteps){
      // Rhs storage
      cds->rhsS    = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
      cds->Sd      = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
      cds->resS    = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
      cds->rhsSd   = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));        
      // Set coefficients equal to flow side time stepper
      cds->SNrk    = cds->fSolver->SNrk;
      cds->Srka    = cds->fSolver->Srka; 
      cds->Srkb    = cds->fSolver->Srkb; 
      cds->Srkc    = cds->fSolver->Srkc;  
      cds->Ue      = cds->fSolver->Ue;
    }


    // MEMORY ALLOCATION
    cds->o_U     = cds->fSolver->o_U; // Set Velocity Field equal to Flow Side
    cds->o_S     = mesh->device.malloc(cds->NSfields*cds->Nstages*Ntotal*sizeof(dfloat), cds->S);
    cds->o_rhsS  = mesh->device.malloc(cds->NSfields*                 Ntotal*sizeof(dfloat), cds->rhsS);
    cds->o_NS    = mesh->device.malloc(cds->NSfields*(cds->Nstages+1)*Ntotal*sizeof(dfloat), cds->NS);
    cds->o_rkS   = mesh->device.malloc(                               Ntotal*sizeof(dfloat), cds->rkS);


    options.getArgs("SBAR", cds->sbar);
    // Copy Kernel Info   
    occa::properties kernelInfo = fSolver->kernelInfo;

    kernelInfo["defines/" "p_alf"]      = cds->alf;
    kernelInfo["defines/" "p_sbar"]     = cds->sbar;
    kernelInfo["defines/" "p_NSfields"] = cds->NSfields;
    kernelInfo["defines/" "p_NVfields"] = cds->NVfields;
    
    for (int r=0;r<mesh->size;r++) {
      if (r==mesh->rank) {
	if (cds->dim==2){ 
	  cds->setScalarFieldKernel =  mesh->device.buildKernel(DCDS "/okl/cdsSetScalarField2D.okl", "cdsSetScalarField2D", kernelInfo);
	}else{
	  cds->setScalarFieldKernel =  mesh->device.buildKernel(DCDS "/okl/cdsSetScalarField3D.okl", "cdsSetScalarField3D", kernelInfo);  
	}
      }
      MPI_Barrier(mesh->comm);
    }

    // Set scalar field
    cds->setScalarFieldKernel(mesh->Nelements,
			      cds->startTime,
			      mesh->o_x,
			      mesh->o_y,
			      mesh->o_z,
			      cds->sOffset,
			      cds->o_S);
    cds->o_S.copyTo(cds->S);
    cds->options = options; 
  cds->options.setArgs("KRYLOV SOLVER",        options.getArgs("SCALAR KRYLOV SOLVER"));
  cds->options.setArgs("DISCRETIZATION",       options.getArgs("SCALAR DISCRETIZATION"));
  cds->options.setArgs("BASIS",                options.getArgs("SCALAR BASIS"));
  cds->options.setArgs("PRECONDITIONER",       options.getArgs("SCALAR PRECONDITIONER"));
  cds->options.setArgs("MULTIGRID COARSENING", options.getArgs("SCALAR MULTIGRID COARSENING"));
  cds->options.setArgs("MULTIGRID SMOOTHER",   options.getArgs("SCALAR MULTIGRID SMOOTHER"));
  cds->options.setArgs("PARALMOND CYCLE",      options.getArgs("SCALAR PARALMOND CYCLE"));
  cds->options.setArgs("PARALMOND SMOOTHER",   options.getArgs("SCALAR PARALMOND SMOOTHER"));
  cds->options.setArgs("PARALMOND PARTITION",  options.getArgs("SCALAR PARALMOND PARTITION"));

  if (mesh->rank==0) printf("==================ELLIPTIC SOLVE SETUP=========================\n");

  // SetUp Boundary Flags types for Elliptic Solve
  // bc = 1 -> wall bc = 2 -> inflow bc = 3 -> outflow
  // bc = 4 -> x-aligned slip, bc = 5 -> y-aligned slip, bc = 6 -> z-aligned slip

  int sBCType[7] = {0,1,1,2,1,1,1}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
 
  //Solver tolerances 
  cds->TOL = 1E-8;

  // Use third Order Velocity Solve: full rank should converge for low orders
  if (mesh->rank==0) printf("==================VELOCITY SOLVE SETUP=========================\n");

  cds->solver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  cds->solver->mesh = mesh;
  cds->solver->options = cds->options;
  cds->solver->dim = cds->dim;
  cds->solver->elementType = cds->elementType;
  cds->solver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(cds->solver->BCType,sBCType,7*sizeof(int));
  ellipticSolveSetup(cds->solver, cds->lambda, kernelInfo); 

  //make node-wise boundary flags
  cds->mapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) cds->mapB[n+e*mesh->Np] = 1E9;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          cds->mapB[fid+e*mesh->Np] = mymin(bc,cds->mapB[fid+e*mesh->Np]);
        }
      }
    }
  }

  ogsGatherScatter(cds->mapB, ogsInt, ogsMin, mesh->ogs);
  
  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (cds->mapB[n] == 1E9) {
      cds->mapB[n] = 0.;
    }
  }
  cds->o_mapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), cds->mapB);

  kernelInfo["defines/" "p_blockSize"]= blockSize;
  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  if(options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    kernelInfo["defines/" "p_EXTBDF"]= 1;
  else
    kernelInfo["defines/" "p_EXTBDF"]= 0;

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = mymax(1,256/mesh->Np); // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1,256/maxNodes); // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;


  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo["defines/" "p_maxNodesVolumeCub"]= maxNodesVolumeCub;
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo["defines/" "p_maxNodesSurfaceCub"]=maxNodesSurfaceCub;
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1); // works for CUDA
  //
  kernelInfo["defines/" "p_cubNblockV"]=cubNblockV;
  kernelInfo["defines/" "p_cubNblockS"]=cubNblockS;

  


#if 1
  if (mesh->rank==0) {
    printf("maxNodes: %d \t  NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
    printf("maxNodesVolCub: %d \t maxNodesSurCub: %d \t NblockVCub: %d \t NblockSCub: %d  \n", maxNodesVolumeCub,maxNodesSurfaceCub, cubNblockV, cubNblockS);

    printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);
  }
#endif
  

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};

    cds->o_rkC  = mesh->device.malloc(4*sizeof(dfloat),rkC);
    cds->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
    cds->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
    cds->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 

    cds->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 

    // cds->o_prkA = cds->o_extbdfC;
    // cds->o_prkB = cds->o_extbdfC;
  }else{
    printf("Only BDF is implemented\n");
    exit(EXIT_FAILURE);
  }

  // MEMORY ALLOCATION
  cds->o_rhsS  = mesh->device.malloc(cds->NSfields*                 Ntotal*sizeof(dfloat), cds->rhsS);
  cds->o_NS    = mesh->device.malloc(cds->NSfields*(cds->Nstages+1)*Ntotal*sizeof(dfloat), cds->NS);
  cds->o_rkS   = mesh->device.malloc(              Ntotal*sizeof(dfloat), cds->rkS);  

  if(mesh->totalHaloPairs){//halo setup
    dlong haloBytes = mesh->totalHaloPairs*mesh->Np*(cds->NSfields)*sizeof(dfloat);
    dlong gatherBytes = cds->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat);
    cds->o_haloBuffer = mesh->device.malloc(haloBytes);

#if 0
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(haloBytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(haloBytes, NULL);
    occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc(gatherBytes, NULL);
    
    cds->sendBuffer    = (dfloat*) o_sendBuffer.getMappedPointer();
    cds->recvBuffer    = (dfloat*) o_recvBuffer.getMappedPointer();
    cds->haloGatherTmp = (dfloat*) o_gatherTmpPinned.getMappedPointer();
#endif
    occa::memory o_sendBuffer, o_recvBuffer,o_gatherTmpPinned;

    cds->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, cds->o_sendBuffer);
    cds->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, cds->o_recvBuffer);

    cds->haloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, gatherBytes, NULL, cds->o_gatherTmpPinned);
    
    cds->o_haloGatherTmp = mesh->device.malloc(gatherBytes,  cds->haloGatherTmp);
  }

  // set kernel name suffix
  char *suffix;
  
  if(cds->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(cds->elementType==QUADRILATERALS){
    if(cds->dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D"); 
  }
  if(cds->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(cds->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r=0;r<mesh->size;r++) {
    if (r==mesh->rank) {
      
      sprintf(fileName, DCDS "/okl/cdsHaloExchange.okl");
      sprintf(kernelName, "cdsHaloExtract");
      cds->haloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      sprintf(kernelName, "cdsHaloScatter");
      cds->haloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      // Probably we dont need this kernel !!!!!!
      // if(cds->dim==3 && cds->elementType==QUADRILATERALS){
      //	sprintf(fileName, DCDS "/okl/insConstrainQuad3D.okl");
      //	sprintf(kernelName, "insConstrainQuad3D");
      //	cds->constrainKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      // }
      
      // ===========================================================================
      
      sprintf(fileName, DCDS "/okl/cdsAdvection%s.okl", suffix);

      // needed to be implemented
      sprintf(kernelName, "cdsAdvectionCubatureVolume%s", suffix);
      cds->advectionCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

     
      sprintf(kernelName, "cdsAdvectionCubatureSurface%s", suffix);
      cds->advectionCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      sprintf(kernelName, "cdsAdvectionVolume%s", suffix);
      cds->advectionVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsAdvectionSurface%s", suffix);
      cds->advectionSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);


      // ===========================================================================
      
      sprintf(fileName, DCDS "/okl/cdsHelmholtzRhs%s.okl", suffix);
      sprintf(kernelName, "cdsHelmholtzRhsEXTBDF%s", suffix);
      cds->helmholtzRhsKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      sprintf(fileName, DCDS "/okl/cdsHelmholtzBC%s.okl", suffix);
      sprintf(kernelName, "cdsHelmholtzIpdgBC%s", suffix);
      cds->helmholtzRhsIpdgBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsHelmholtzBC%s", suffix);
      cds->helmholtzRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsHelmholtzAddBC%s", suffix);
      cds->helmholtzAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      /*
      // ===========================================================================
      if(cds->dim==3 && cds->options.compareArgs("OUTPUT TYPE","ISO")){
      sprintf(fileName, DCDS "/okl/insIsoSurface3D.okl");
      sprintf(kernelName, "insIsoSurface3D");

      cds->isoSurfaceKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  
      }
      */

      // Not implemented for Quad 3D yet !!!!!!!!!!
      if(cds->Nsubsteps){
      // Note that resU and resV can be replaced with already introduced buffer
      cds->o_Ue     = mesh->device.malloc(cds->NVfields*Ntotal*sizeof(dfloat), cds->Ue);
      cds->o_Sd     = mesh->device.malloc(cds->NSfields*Ntotal*sizeof(dfloat), cds->Sd);
      cds->o_resS   = mesh->device.malloc(cds->NSfields*Ntotal*sizeof(dfloat), cds->resS);
      cds->o_rhsSd  = mesh->device.malloc(cds->NSfields*Ntotal*sizeof(dfloat), cds->rhsSd);

      //if(cds->elementType==HEXAHEDRA)
      // cds->o_cSd = mesh->device.malloc(cds->NSfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), cds->cSd);
      //else 
      // cds->o_cSd = cds->o_Sd;

      sprintf(fileName, DHOLMES "/okl/scaledAdd.okl");
      sprintf(kernelName, "scaledAddwOffset");
      cds->scaledAddKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

       sprintf(fileName, DCDS "/okl/cdsSubCycle%s.okl", suffix);
       sprintf(kernelName, "cdsSubCycleVolume%s", suffix);
       cds->subCycleVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsSubCycleSurface%s", suffix);
      cds->subCycleSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

     sprintf(kernelName, "cdsSubCycleCubatureVolume%s", suffix);
     cds->subCycleCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsSubCycleCubatureSurface%s", suffix);
      cds->subCycleCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(fileName, DCDS "/okl/cdsSubCycle.okl");
      sprintf(kernelName, "cdsSubCycleRKUpdate");
      cds->subCycleRKUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "cdsSubCycleExt");
      cds->subCycleExtKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }
   
    }
    MPI_Barrier(mesh->comm);
  }
#endif      
}
 #if 0
  //sanity checking
  if (options.compareArgs("ELEMENT TYPE","HEXAHEDRA"){
      printf("ERROR: Only HEXAHEDRAL element is available currently\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  
    dlong Ntotal = mesh->Np*mesh->Nelements;
    dlong Nblock = mymax(1,(Ntotal+blockSize-1)/blockSize);
    dlong Nhalo  = mesh->Np*mesh->totalHaloPairs;
    dlong Nall   = Ntotal + Nhalo;

    // Assumes NStages are initiliazed in solver side.... 
    cds->S     = (dfloat*) calloc(cds->NSfields*cds->Nstages*Ntotal,sizeof(dfloat));
    //additional field storage
    cds->NS    = (dfloat*) calloc(cds->NSfields*(cds->Nstages+1)*Ntotal,sizeof(dfloat));
    cds->rkS   = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
    cds->rkNS  = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));

    if(cds->Nsubsteps){
      // Rhs storage
      cds->rhsS    = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
      cds->Sd      = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
      cds->resS    = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
      cds->rhsSd   = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));        
      // Set coefficients equal to flow side
      cds->SNrk    = cds->fSolver->SNrk;
      cds->Srka    = cds->fSolver->Srka; 
      cds->Srkb    = cds->fSolver->Srkb; 
      cds->Srkc    = cds->fSolver->Srkc;  
      cds->Ue      = cds->fSolver->Ue;
    }
    
    // MEMORY ALLOCATION
    cds->o_rhsS  = mesh->device.malloc(cds->NSfields*                 Ntotal*sizeof(dfloat), cds->rhsS);
    cds->o_NS    = mesh->device.malloc(cds->NSfields*(cds->Nstages+1)*Ntotal*sizeof(dfloat), cds->NS);
    cds->o_rkS   = mesh->device.malloc(              Ntotal*sizeof(dfloat), cds->rkS);  
    cds->o_rkNS  = mesh->device.malloc(cds->NVfields*Ntotal*sizeof(dfloat), cds->rkNS);

    cds->rho = cds->fSolver->rho;
    cds->nu  = cds->fSolver->nu; 
    //Reynolds number
    cds->Re     = cds->fSolver->Re;          // Reynolds number
    cds->Pr     = cds->nu/cds->alf;          // Prandtl number

    // Get Solver parameters
    options.getArgs("CONDUCTIVITY", cds->k);
    options.getArgs("SPECIFIC HEAT", cds->cp);
    
    cds->alf    = cds->k/(cds->rho*cds->cp); // Thermal diffusivity
    cds->ialf   = 1.0/cds->alf;              // Inverse diff. 
    
    cds->o_U    = cds->fSolver->o_U; 
    cds->o_S    = mesh->device.malloc(cds->NSfields*cds->Nstages*Ntotal*sizeof(dfloat), cds->S);

    for (int r=0;r<mesh->size;r++) {
      if (r==mesh->rank) {
	if (cds->dim==2){ 
	  cds->setScalarFieldKernel =  mesh->device.buildKernel(DCDS "/okl/cdsSetScalarField2D.okl", "cdsSetScalarField2D", kernelInfo);
	}else{
	  cds->setScalarFieldKernel =  mesh->device.buildKernel(DCDS "/okl/cdsSetScalarField3D.okl", "cdsSetScalarField3D", kernelInfo);  
	}
      }
      MPI_Barrier(mesh->comm);
    }

    // Set scalar field
    cds->setScalarFieldKernel(mesh->Nelements,
			      cds->startTime,
			      mesh->o_x,
			      mesh->o_y,
			      mesh->o_z,
			      cds->sOffset,
			      cds->o_S);
    cds->o_S.copyTo(cds->S);

    // save initial time-step estimate 
    cds->dt = cds->fSolver->dt; 
    cds->g0 = cds->fSolver->g0; 
    cds->idt     = 1.0/cds->dt;
    cds->lambda  = cds->g0 / (cds->dt * cds->alf);


    kernelInfo["defines/" "p_alf"]     = cds->alf;


    // Set time integrator the same with flow solver
    if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
      cds->o_rkC     = cds->fSolver->o_rkC; 
      cds->o_extbdfA = cds->fSolver->o_extbdfA; 
      cds->o_extbdfB = cds->fSolver->o_extbdfB; 
      cds->o_extbdfC = cds->fSolver->o_extbdfC; 
      cds->o_extC    = cds->fSolver->o_extC; 
    }else{
      printf("Only BDF is implemented\n");
      exit(EXIT_FAILURE);
    }


    if(mesh->totalHaloPairs){//halo setup
      dlong haloBytes = mesh->totalHaloPairs*mesh->Np*(cds->NSfields)*sizeof(dfloat);
      dlong gatherBytes = cds->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat);
      cds->o_haloBuffer = mesh->device.malloc(haloBytes);

#if 0
      occa::memory o_sendBuffer = mesh->device.mappedAlloc(haloBytes, NULL);
      occa::memory o_recvBuffer = mesh->device.mappedAlloc(haloBytes, NULL);
      occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc(gatherBytes, NULL);
    
      cds->sendBuffer    = (dfloat*) o_sendBuffer.getMappedPointer();
      cds->recvBuffer    = (dfloat*) o_recvBuffer.getMappedPointer();
      cds->haloGatherTmp = (dfloat*) o_gatherTmpPinned.getMappedPointer();
#endif
      occa::memory o_sendBuffer, o_recvBuffer,o_gatherTmpPinned;

      cds->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, cds->o_sendBuffer);
      cds->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, haloBytes, NULL, cds->o_recvBuffer);

      cds->haloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, gatherBytes, NULL, cds->o_gatherTmpPinned);
    
      cds->o_haloGatherTmp = mesh->device.malloc(gatherBytes,  cds->haloGatherTmp);
    }











    
    

    //setup async halo stream
    cds->defaultStream = mesh->defaultStream;
    cds->dataStream    = mesh->dataStream;

    dlong Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    if(Nbytes>0){
#if 0
      occa::memory o_sendBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
      occa::memory o_recvBuffer = mesh->device.mappedAlloc(Nbytes, NULL);

      cds->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
      cds->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();


      occa::memory o_gradSendBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);
      occa::memory o_gradRecvBuffer = mesh->device.mappedAlloc(2*Nbytes, NULL);

      cds->gradSendBuffer = (dfloat*) o_gradSendBuffer.getMappedPointer();
      cds->gradRecvBuffer = (dfloat*) o_gradRecvBuffer.getMappedPointer();
#endif

      cds->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, Nbytes, NULL, cds->o_sendBuffer);
      cds->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, Nbytes, NULL, cds->o_recvBuffer);
      cds->gradSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, 2*Nbytes, NULL, cds->o_gradSendBuffer);
      cds->gradRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, 2*Nbytes, NULL, cds->o_gradRecvBuffer);

    }else{
      cds->sendBuffer = NULL;
      cds->recvBuffer = NULL;
    }
    mesh->device.setStream(cds->defaultStream);

    cds->type = strdup(dfloatString);

    cds->Nblock = Nblock;
    cds->Nblock2 = Nblock2;

    //fill geometric factors in halo
    if(mesh->totalHaloPairs){
      dlong Nlocal = mesh->Nelements;
      dlong Nhalo = mesh->totalHaloPairs;
      size_t Nbytes = mesh->Nvgeo*sizeof(dfloat);

      if (cds->elementType==QUADRILATERALS || cds->elementType==HEXAHEDRA) {
	Nlocal *= mesh->Np;
	Nhalo *= mesh->Np;
	Nbytes *= mesh->Np;
      }

      dfloat *vgeoSendBuffer = (dfloat*) calloc(Nhalo*mesh->Nvgeo, sizeof(dfloat));

      // import geometric factors from halo elements
      mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat));

      meshHaloExchange(mesh,
		       Nbytes,
		       mesh->vgeo,
		       vgeoSendBuffer,
		       mesh->vgeo + Nlocal*mesh->Nvgeo);

      mesh->o_vgeo =
	mesh->device.malloc((Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);
      free(vgeoSendBuffer);
    }


    //build inverse of mass matrix
    mesh->invMM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    for (int n=0;n<mesh->Np*mesh->Np;n++)
      mesh->invMM[n] = mesh->MM[n];
    matrixInverse(mesh->Np,mesh->invMM);


    //check all the bounaries for a Dirichlet
    bool allNeumann = (lambda==0) ? true :false;
    cds->allNeumannPenalty = 1.;
    hlong localElements = (hlong) mesh->Nelements;
    hlong totalElements = 0;
    MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, mesh->comm);
    cds->allNeumannScale = 1./sqrt((dfloat)mesh->Np*totalElements);

    cds->EToB = (int *) calloc(mesh->Nelements*mesh->Nfaces,sizeof(int));
    for (dlong e=0;e<mesh->Nelements;e++) {
      for (int f=0;f<mesh->Nfaces;f++) {
	int bc = mesh->EToB[e*mesh->Nfaces+f];
	if (bc>0) {
	  int BC = cds->BCType[bc]; //get the type of boundary
	  cds->EToB[e*mesh->Nfaces+f] = BC; //record it
	  if (BC!=2) allNeumann = false; //check if its a Dirchlet
	}
      }
    }


    // !!!!!! Removed MPI::BOOL since some mpi versions complains about it !!!!!
    int lallNeumann, gallNeumann;
    lallNeumann = allNeumann ? 0:1;
    MPI_Allreduce(&lallNeumann, &gallNeumann, 1, MPI_INT, MPI_SUM, mesh->comm);
    cds->allNeumann = (gallNeumann>0) ? false: true;

    // MPI_Allreduce(&allNeumann, &(cds->allNeumann), 1, MPI::BOOL, MPI_LAND, mesh->comm);
    if (mesh->rank==0&& options.compareArgs("VERBOSE","TRUE")) printf("allNeumann = %d \n", cds->allNeumann);

    //set surface mass matrix for continuous boundary conditions
    mesh->sMT = (dfloat *) calloc(mesh->Np*mesh->Nfaces*mesh->Nfp,sizeof(dfloat));
    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Nfp*mesh->Nfaces;m++) {
	dfloat MSnm = 0;
	for (int i=0;i<mesh->Np;i++){
	  MSnm += mesh->MM[n+i*mesh->Np]*mesh->LIFT[m+i*mesh->Nfp*mesh->Nfaces];
	}
	mesh->sMT[n+m*mesh->Np]  = MSnm;
      }
    }
    mesh->o_sMT = mesh->device.malloc(mesh->Np*mesh->Nfaces*mesh->Nfp*sizeof(dfloat), mesh->sMT);

    //copy boundary flags
    cds->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(int), cds->EToB);

#if 0
    if (mesh->rank==0 && options.compareArgs("VERBOSE","TRUE"))
      occa::setVerboseCompilation(true);
    else
      occa::setVerboseCompilation(false);
#endif

    //setup an unmasked gs handle
    int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;
    meshParallelGatherScatterSetup(mesh, Ntotal, mesh->globalIds, mesh->comm, verbose);

    //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
    cds->mapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
    for (dlong e=0;e<mesh->Nelements;e++) {
      for (int n=0;n<mesh->Np;n++) cds->mapB[n+e*mesh->Np] = 1E9;
      for (int f=0;f<mesh->Nfaces;f++) {
	int bc = mesh->EToB[f+e*mesh->Nfaces];
	if (bc>0) {
	  for (int n=0;n<mesh->Nfp;n++) {
	    int BCFlag = cds->BCType[bc];
	    int fid = mesh->faceNodes[n+f*mesh->Nfp];
	    cds->mapB[fid+e*mesh->Np] = mymin(BCFlag,cds->mapB[fid+e*mesh->Np]);
	  }
	}
      }
    }
    ogsGatherScatter(cds->mapB, ogsInt, ogsMin, mesh->ogs);

    //use the bc flags to find masked ids
    cds->Nmasked = 0;
    for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
      if (cds->mapB[n] == 1E9) {
	cds->mapB[n] = 0.;
      } else if (cds->mapB[n] == 1) { //Dirichlet boundary
	cds->Nmasked++;
      }
    }

  
    cds->o_mapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), cds->mapB);

    cds->maskIds = (dlong *) calloc(cds->Nmasked, sizeof(dlong));
    cds->Nmasked =0; //reset
    for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
      if (cds->mapB[n] == 1) cds->maskIds[cds->Nmasked++] = n;
    }
    if (cds->Nmasked) cds->o_maskIds = mesh->device.malloc(cds->Nmasked*sizeof(dlong), cds->maskIds);

    //make a masked version of the global id numbering
    mesh->maskedGlobalIds = (hlong *) calloc(Ntotal,sizeof(hlong));
    memcpy(mesh->maskedGlobalIds, mesh->globalIds, Ntotal*sizeof(hlong));
    for (dlong n=0;n<cds->Nmasked;n++)
      mesh->maskedGlobalIds[cds->maskIds[n]] = 0;

    //use the masked ids to make another gs handle
    cds->ogs = ogsSetup(Ntotal, mesh->maskedGlobalIds, mesh->comm, verbose, mesh->device);
    cds->o_invDegree = cds->ogs->o_invDegree;

    /*preconditioner setup */
    cds->precon = (precon_t*) calloc(1, sizeof(precon_t));

    kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

    if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
      kernelInfo["compiler_flags"] += "-Xptxas -dlcm=ca";
    }

    if(mesh->device.mode()=="Serial")
      kernelInfo["compiler_flags"] += "-g";

    // set kernel name suffix
    char *suffix;

    if(cds->elementType==TRIANGLES){
      if(cds->dim==2)
	suffix = strdup("Tri2D");
      else
	suffix = strdup("Tri3D");
    }
    if(cds->elementType==QUADRILATERALS){
      if(cds->dim==2)
	suffix = strdup("Quad2D");
      else
	suffix = strdup("Quad3D");
    }
    if(cds->elementType==TETRAHEDRA)
      suffix = strdup("Tet3D");
    if(cds->elementType==HEXAHEDRA)
      suffix = strdup("Hex3D");

    char fileName[BUFSIZ], kernelName[BUFSIZ];


    for (int r=0;r<mesh->size;r++) {
      if (r==mesh->rank) {

	//mesh kernels
	mesh->haloExtractKernel =
	  mesh->device.buildKernel(DHOLMES "/okl/meshHaloExtract2D.okl",
				   "meshHaloExtract2D",
				   kernelInfo);

	mesh->addScalarKernel =
	  mesh->device.buildKernel(DHOLMES "/okl/addScalar.okl",
				   "addScalar",
				   kernelInfo);

	mesh->maskKernel =
	  mesh->device.buildKernel(DHOLMES "/okl/mask.okl",
				   "mask",
				   kernelInfo);


	kernelInfo["defines/" "p_blockSize"]= blockSize;


	mesh->sumKernel =
	  mesh->device.buildKernel(DHOLMES "/okl/sum.okl",
				   "sum",
				   kernelInfo);

	cds->weightedInnerProduct1Kernel =
	  mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct1.okl",
				   "weightedInnerProduct1",
				   kernelInfo);

	cds->weightedInnerProduct2Kernel =
	  mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl",
				   "weightedInnerProduct2",
				   kernelInfo);

	cds->innerProductKernel =
	  mesh->device.buildKernel(DHOLMES "/okl/innerProduct.okl",
				   "innerProduct",
				   kernelInfo);

	cds->weightedNorm2Kernel =
	  mesh->device.buildKernel(DHOLMES "/okl/weightedNorm2.okl",
				   "weightedNorm2",
				   kernelInfo);

	cds->norm2Kernel =
	  mesh->device.buildKernel(DHOLMES "/okl/norm2.okl",
				   "norm2",
				   kernelInfo);


	cds->scaledAddKernel =
          mesh->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
				   "scaledAdd",
				   kernelInfo);

	cds->dotMultiplyKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
				   "dotMultiply",
				   kernelInfo);

	cds->dotDivideKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotDivide.okl",
				   "dotDivide",
				   kernelInfo);

	// add custom defines
	kernelInfo["defines/" "p_NpP"]= (mesh->Np+mesh->Nfp*mesh->Nfaces);
	kernelInfo["defines/" "p_Nverts"]= mesh->Nverts;

	//sizes for the coarsen and prolongation kernels. degree N to degree 1
	kernelInfo["defines/" "p_NpFine"]= mesh->Np;
	kernelInfo["defines/" "p_NpCoarse"]= mesh->Nverts;


	if (cds->elementType==QUADRILATERALS || cds->elementType==HEXAHEDRA) {
	  kernelInfo["defines/" "p_NqFine"]= mesh->N+1;
	  kernelInfo["defines/" "p_NqCoarse"]= 2;
	}

	kernelInfo["defines/" "p_NpFEM"]= mesh->NpFEM;

	int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
	kernelInfo["defines/" "p_Nmax"]= Nmax;

	int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
	kernelInfo["defines/" "p_maxNodes"]= maxNodes;

	int NblockV = mymax(1,maxNthreads/mesh->Np); // works for CUDA
	int NnodesV = 1; //hard coded for now
	kernelInfo["defines/" "p_NblockV"]= NblockV;
	kernelInfo["defines/" "p_NnodesV"]= NnodesV;
	kernelInfo["defines/" "p_NblockVFine"]= NblockV;
	kernelInfo["defines/" "p_NblockVCoarse"]= NblockV;

	int NblockS = mymax(1,maxNthreads/maxNodes); // works for CUDA
	kernelInfo["defines/" "p_NblockS"]= NblockS;

	int NblockP = mymax(1,maxNthreads/(4*mesh->Np)); // get close to maxNthreads threads
	kernelInfo["defines/" "p_NblockP"]= NblockP;

	int NblockG;
	if(mesh->Np<=32) NblockG = ( 32/mesh->Np );
	else NblockG = maxNthreads/mesh->Np;
	kernelInfo["defines/" "p_NblockG"]= NblockG;

	//add standard boundary functions
	char *boundaryHeaderFileName;
	if (cds->dim==2)
	  boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary2D.h");
	else if (cds->dim==3)
	  boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary3D.h");
	kernelInfo["includes"] += boundaryHeaderFileName;


	sprintf(fileName,  DELLIPTIC "/okl/ellipticAx%s.okl", suffix);
	sprintf(kernelName, "ellipticAx%s", suffix);

	occa::properties dfloatKernelInfo = kernelInfo;
	occa::properties floatKernelInfo = kernelInfo;
	floatKernelInfo["defines/" "pfloat"]= "float";
	dfloatKernelInfo["defines/" "pfloat"]= dfloatString;

	cds->AxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);

	if(cds->elementType!=HEXAHEDRA){
	  sprintf(kernelName, "ellipticPartialAx%s", suffix);
	}
	else{
	  if(cds->options.compareArgs("ELEMENT MAP", "TRILINEAR")){
	    sprintf(kernelName, "ellipticPartialAxTrilinear%s", suffix);
	  }else{
	    sprintf(kernelName, "ellipticPartialAx%s", suffix);
	  }
	}

	cds->partialAxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);

	cds->partialFloatAxKernel = mesh->device.buildKernel(fileName,kernelName,floatKernelInfo);

	// Not implemented for Quad3D !!!!!
	if (options.compareArgs("BASIS","BERN")) {

	  sprintf(fileName, DELLIPTIC "/okl/ellipticGradientBB%s.okl", suffix);
	  sprintf(kernelName, "ellipticGradientBB%s", suffix);

	  cds->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	  sprintf(kernelName, "ellipticPartialGradientBB%s", suffix);
	  cds->partialGradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	  sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdgBB%s.okl", suffix);
	  sprintf(kernelName, "ellipticAxIpdgBB%s", suffix);
	  cds->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	  sprintf(kernelName, "ellipticPartialAxIpdgBB%s", suffix);
	  cds->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	} else if (options.compareArgs("BASIS","NODAL")) {

	  sprintf(fileName, DELLIPTIC "/okl/ellipticGradient%s.okl", suffix);
	  sprintf(kernelName, "ellipticGradient%s", suffix);

	  cds->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	  sprintf(kernelName, "ellipticPartialGradient%s", suffix);
	  cds->partialGradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	  sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdg%s.okl", suffix);
	  sprintf(kernelName, "ellipticAxIpdg%s", suffix);
	  cds->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	  sprintf(kernelName, "ellipticPartialAxIpdg%s", suffix);
	  cds->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
	}

	// Use the same kernel with quads for the following kenels
	if(cds->dim==3){
	  if(cds->elementType==QUADRILATERALS)
	    suffix = strdup("Quad2D");
	  else if(cds->elementType==TRIANGLES)
	    suffix = strdup("Tri2D");
	}
	
	sprintf(fileName, DELLIPTIC "/okl/ellipticPreconCoarsen%s.okl", suffix);
	sprintf(kernelName, "ellipticPreconCoarsen%s", suffix);
	cds->precon->coarsenKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	sprintf(fileName, DELLIPTIC "/okl/ellipticPreconProlongate%s.okl", suffix);
	sprintf(kernelName, "ellipticPreconProlongate%s", suffix);
	cds->precon->prolongateKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      

	sprintf(fileName, DELLIPTIC "/okl/ellipticBlockJacobiPrecon.okl");
	sprintf(kernelName, "ellipticBlockJacobiPrecon");
	cds->precon->blockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	sprintf(kernelName, "ellipticPartialBlockJacobiPrecon");
	cds->precon->partialblockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	sprintf(fileName, DELLIPTIC "/okl/ellipticPatchSolver.okl");
	sprintf(kernelName, "ellipticApproxBlockJacobiSolver");
	cds->precon->approxBlockJacobiSolverKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

	if (   cds->elementType == TRIANGLES
	       || cds->elementType == TETRAHEDRA) {
	  cds->precon->SEMFEMInterpKernel =
	    mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMInterp.okl",
				     "ellipticSEMFEMInterp",
				     kernelInfo);

	  cds->precon->SEMFEMAnterpKernel =
	    mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMAnterp.okl",
				     "ellipticSEMFEMAnterp",
				     kernelInfo);
	}
      }
      MPI_Barrier(mesh->comm);
    }

    long long int pre = mesh->device.memoryAllocated();

    occaTimerTic(mesh->device,"PreconditionerSetup");
    ellipticPreconditionerSetup(elliptic, cds->ogs, lambda);
    occaTimerToc(mesh->device,"PreconditionerSetup");

    long long int usedBytes = mesh->device.memoryAllocated()-pre;

    cds->precon->preconBytes = usedBytes;
    }

}
#endif
