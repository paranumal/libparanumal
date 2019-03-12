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

#include "ins.h"
#include "omp.h"
#include <unistd.h>

void insSetScalarSolver(ins_t *ins, setupAide options, occa::properties &kernelInfo){

// int rank, size; 
// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// MPI_Comm_size(MPI_COMM_WORLD, &size); 

mesh_t *mesh = ins->mesh; 
ins->sSolver = (cds_t*) calloc(1, sizeof(cds_t));

cds_t *cds = ins->sSolver; 

if(mesh->rank==0) printf("......Setting Scalar Solver........ \n");

if(mesh->rank==0) printf("Setting Primitives.....\n");
// set mesh, options
cds->mesh        = mesh; 
cds->options     = options; 
cds->elementType = ins->elementType; 
cds->dim         = ins->dim; 
cds->NVfields    = ins->NVfields;
// Number of scalar field is hard coded 
cds->NSfields    = 1; 

if(mesh->rank==0) printf("Setting Time Stepper Info.....\n");
cds->extbdfA = ins->extbdfA;
cds->extbdfB = ins->extbdfB;
cds->extbdfC = ins->extbdfC;
cds->extC    = ins->extC   ;
//cds->extC    = (dfloat*) calloc(3, sizeof(dfloat));
//
cds->Nstages       = ins->Nstages; 
cds->temporalOrder = ins->temporalOrder; 
cds->g0            = ins->g0; 


dlong Nlocal = mesh->Np*mesh->Nelements;
dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

cds->Ntotal      = Ntotal;
cds->vOffset     = Ntotal; // keep it different for now.....b
cds->sOffset     = Ntotal;
cds->Nblock      = (Nlocal+blockSize-1)/blockSize;

// Solution storage at interpolation nodes
cds->U     = (dfloat*) calloc(cds->NVfields*(cds->Nstages+0)*Ntotal,sizeof(dfloat));
cds->S     = (dfloat*) calloc(cds->NSfields*(cds->Nstages+0)*Ntotal,sizeof(dfloat));
cds->NS    = (dfloat*) calloc(cds->NSfields*(cds->Nstages+1)*Ntotal,sizeof(dfloat));
cds->rkS   = (dfloat*) calloc(cds->NSfields                 *Ntotal,sizeof(dfloat));


cds->U = ins->U; // Point to INS side Velocity

cds->S     = (dfloat*) calloc(cds->NSfields*(cds->Nstages+0)*Ntotal,sizeof(dfloat));
cds->NS    = (dfloat*) calloc(cds->NSfields*(cds->Nstages+1)*Ntotal,sizeof(dfloat));
cds->rkS   = (dfloat*) calloc(cds->NSfields                 *Ntotal,sizeof(dfloat));
cds->rhsS  = (dfloat*) calloc(cds->NSfields                 *Ntotal,sizeof(dfloat));
// Use Nsubsteps if INS does to prevent stability issues
cds->Nsubsteps = ins->Nsubsteps; 


if(cds->Nsubsteps){
	// This memory can be reduced, check later......!!!!!!!
    cds->Sd      = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
    cds->resS    = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));
    cds->rhsSd   = (dfloat*) calloc(cds->NSfields*Ntotal,sizeof(dfloat));        
    cds->Ue      = (dfloat*) calloc(cds->NVfields*Ntotal,sizeof(dfloat));
    // 
    cds->SNrk    = ins->SNrk;  
    //
    cds->Srka = ins->Srka; 
    cds->Srkb = ins->Srkb; 
    cds->Srkc = ins->Srkc; 
  }

 
 options.getArgs("UBAR", cds->ubar);
  options.getArgs("VBAR", cds->vbar);
  if (cds->dim==3)
    options.getArgs("WBAR", cds->wbar);
  options.getArgs("SBAR", cds->sbar);

  options.getArgs("CONDUCTIVITY", cds->k);
  options.getArgs("SPECIFIC HEAT", cds->cp);
  options.getArgs("DENSITY", cds->rho);
  options.getArgs("VISCOSITY", cds->nu);

  cds->alf    = cds->k/(cds->rho*cds->cp); // Thermal diffusivity
  cds->ialf   = 1.0/cds->alf;              // Inverse diff. 


  //Reynolds number
  cds->Re     = cds->ubar/cds->nu;         // Reynolds number
  cds->Pr     = cds->nu/cds->alf;          // Prandtl number
   

  occa::properties kernelInfoH     = kernelInfo; 
   // ADD-DEFINES
  kernelInfo["defines/" "p_sbar"]= cds->sbar;
  kernelInfo["defines/" "p_ubar"]= cds->ubar;
  kernelInfo["defines/" "p_vbar"]= cds->vbar;
  kernelInfo["defines/" "p_wbar"]= cds->wbar;

  kernelInfo["defines/" "p_NSfields"]= cds->NSfields;
  kernelInfo["defines/" "p_NVfields"]= cds->NVfields;
  kernelInfo["defines/" "p_NTfields"]= (cds->NVfields+cds->NSfields);
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  cds->Nstages;
  kernelInfo["defines/" "p_SUBCYCLING"]=  cds->Nsubsteps;


   //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();

  cds->o_U = ins->o_U; // point to INS velocity very important !!!!!
  cds->o_S = mesh->device.malloc(cds->NSfields*(cds->Nstages+0)*Ntotal*sizeof(dfloat), cds->S);

  printf("Compiling CDS kernels.....\n");
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

  cds->startTime =ins->startTime;

  // options.getArgs("START TIME", cds->startTime);
  cds->setScalarFieldKernel(mesh->Nelements,
			    cds->startTime,
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    cds->sOffset,
			    cds->o_S);

  // cds->o_S.copyTo(cds->S);

  cds->dt  = ins->dt; 
  cds->sdt = ins->sdt; 
  cds->NtimeSteps = ins->NtimeSteps; 

   if (cds->Nsubsteps && mesh->rank==0) printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", cds->dt, cds->sdt, cds->dt/cds->sdt);

  cds->idt     = 1.0/cds->dt;
  cds->lambda  = cds->g0 / (cds->dt * cds->alf);

  kernelInfo["defines/" "p_alf"]     = cds->alf;

  //make option objects for elliptc solvers
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
  if (mesh->rank==0) printf("==================HELMHOLTZ SETUP SOLVE SETUP=========================\n");

  cds->solver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  cds->solver->mesh = mesh;
  cds->solver->options = cds->options;
  cds->solver->dim = cds->dim;
  cds->solver->elementType = cds->elementType;
  cds->solver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(cds->solver->BCType,sBCType,7*sizeof(int));
  ellipticSolveSetup(cds->solver, cds->lambda, kernelInfoH); 

  //make node-wise boundary flags // It is batter hold a different one, 
  // because wall may not give Dirichlet in temperature solve !!!!!!
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

  //

  dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};

	cds->o_rkC     = ins->o_rkC    ;
	cds->o_extbdfA = ins->o_extbdfA;
	cds->o_extbdfB = ins->o_extbdfB;
	cds->o_extbdfC = ins->o_extbdfC; 

	 cds->o_extC = ins->o_extC;
	 // cds->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 

	cds->o_prkA = ins->o_extbdfC;
	cds->o_prkB = ins->o_extbdfC;

// MEMORY ALLOCATION
  cds->o_rhsS  = mesh->device.malloc(cds->NSfields*                 Ntotal*sizeof(dfloat), cds->rhsS);
  cds->o_NS    = mesh->device.malloc(cds->NSfields*(cds->Nstages+1)*Ntotal*sizeof(dfloat), cds->NS);
  cds->o_rkS   = mesh->device.malloc(                               Ntotal*sizeof(dfloat), cds->rkS);  



   if(mesh->totalHaloPairs){//halo setup
    dlong haloBytes   = mesh->totalHaloPairs*mesh->Np*(cds->NSfields + cds->NVfields)*sizeof(dfloat);
    dlong gatherBytes = (cds->NSfields+cds->NVfields)*mesh->ogs->NhaloGather*sizeof(dfloat);
    cds->o_haloBuffer = mesh->device.malloc(haloBytes);

   

    // dlong vhaloBytes   = mesh->totalHaloPairs*mesh->Np*(cds->NSfields)*sizeof(dfloat);
    // dlong vgatherBytes = cds->NSfields*mesh->ogs->NhaloGather*sizeof(dfloat);
    // cds->o_vhaloBuffer = mesh->device.malloc(vhaloBytes);

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

    // Halo exchange for more efficient subcycling 
    if(cds->Nsubsteps){
      dlong shaloBytes   = mesh->totalHaloPairs*mesh->Np*(cds->NSfields)*sizeof(dfloat);
      dlong sgatherBytes = (cds->NSfields)*mesh->ogs->NhaloGather*sizeof(dfloat);
      cds->o_shaloBuffer = mesh->device.malloc(shaloBytes);

      occa::memory o_ssendBuffer, o_srecvBuffer,o_sgatherTmpPinned;

      cds->ssendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, shaloBytes, NULL, cds->o_ssendBuffer);
      cds->srecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, shaloBytes, NULL, cds->o_srecvBuffer);
      cds->shaloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, sgatherBytes, NULL, cds->o_sgatherTmpPinned);
      cds->o_shaloGatherTmp = mesh->device.malloc(sgatherBytes,  cds->shaloGatherTmp);
   }
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
      
      if(cds->Nsubsteps){
        sprintf(kernelName, "cdsScalarHaloExtract");
        cds->scalarHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
        
        sprintf(kernelName, "cdsScalarHaloScatter");
        cds->scalarHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);    
      } 
      
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

    
      if(cds->Nsubsteps){
        // Note that resU and resV can be replaced with already introduced buffer
        // cds->o_Ue     = cds->Ue;
        cds->o_Ue     = mesh->device.malloc(cds->NVfields*Ntotal*sizeof(dfloat), cds->Ue);
        // Can use previous memories
        cds->o_Sd     = mesh->device.malloc(cds->NSfields*Ntotal*sizeof(dfloat), cds->Sd);
        cds->o_resS   = mesh->device.malloc(cds->NSfields*Ntotal*sizeof(dfloat), cds->resS);
        cds->o_rhsSd  = mesh->device.malloc(cds->NSfields*Ntotal*sizeof(dfloat), cds->rhsSd);

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

 if(mesh->rank==0) printf("INS + CDS Setup is done \n");

 #if 0
  printf(" Solver Parameters......\n");
  printf("alfa\t:\t %.8e \n", cds->alf);
  printf("invalfa\t:\t %.8e \n", cds->ialf);
  printf("k\t:\t %.8e \n", cds->k);
  printf("cp\t:\t %.8e \n", cds->cp);
  printf("rho\t:\t %.8e \n", cds->rho);
  printf("nu\t:\t %.8e \n", cds->nu);
  printf("Re\t:\t %.8e \n", cds->Re);
  printf("Pr\t:\t %.8e \n", cds->Pr);
  printf("dt\t:\t %.8e \n", cds->dt);
  printf("sdt\t:\t %.8e \n", cds->sdt);
  printf("invdt\t:\t %.8e \n", cds->idt);
  printf("ialf\t:\t %.8e \n", cds->ialf);
  printf("Nsubsteps\t:\t %02d \n", cds->Nsubsteps);
  printf("Nstages\t:\t %02d \n", cds->Nstages);
  printf("SNrk\t:\t %02d \n", cds->SNrk); 
  printf("ExplicitOrder\t:\t %02d \n", cds->ExplicitOrder);       

#endif


}