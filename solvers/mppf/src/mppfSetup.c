#include "mppf.h"
#include "omp.h"
#include <unistd.h>

mppf_t *mppfSetup(mesh_t *mesh, setupAide options){
	
	mppf_t *mppf = (mppf_t *) calloc(1,sizeof(mppf_t));

	mppf->mesh    = mesh; 
	mppf->options = options; 

	options.getArgs("MESH DIMENSION", mppf->dim);
	options.getArgs("ELEMENT TYPE", mppf->elementType);

	// Flow Side
	mppf->NVfields = (mppf->dim==3) ? 3:2; //  Total Number of Velocity Fields
	mppf->NTfields = (mppf->dim==3) ? 4:3; // Total Velocity + Pressure

	// Phase Tracking Side
	mesh->Nfields = 1; // set mesh Nfields to 1, this may cause problems in meshOccaSetup 


	if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
		mppf->extbdfA = (dfloat*) calloc(3, sizeof(dfloat));
		mppf->extbdfB = (dfloat*) calloc(3, sizeof(dfloat));
		mppf->extbdfC = (dfloat*) calloc(3, sizeof(dfloat));
		mppf->extC    = (dfloat*) calloc(3, sizeof(dfloat));
	}

	if (       options.compareArgs("TIME INTEGRATOR", "EXTBDF1")) {
		mppf->Nstages = 1;
		mppf->temporalOrder = 1;
		mppf->g0 = 1.0;
	} else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF2")) {
		mppf->Nstages = 2;
		mppf->temporalOrder = 2;
		mppf->g0 = 1.5;
	} else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF3")) {
		mppf->Nstages = 3;
		mppf->temporalOrder = 3;
		mppf->g0 = 11.f/6.f;
	}

	// Initialize fields
	dlong Nlocal = mesh->Np*mesh->Nelements;
	dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

	mppf->Ntotal = Ntotal;
	mppf->fieldOffset = Ntotal;
	mppf->Nblock = (Nlocal+blockSize-1)/blockSize;

	// flow side storage
	mppf->U      	= (dfloat*) calloc(mppf->NVfields*  mppf->Nstages*Ntotal,sizeof(dfloat));
	mppf->P      	= (dfloat*) calloc(                 mppf->Nstages*Ntotal,sizeof(dfloat));

	// interface side storage
	mppf->Phi    	= (dfloat*) calloc(                 mppf->Nstages*Ntotal,sizeof(dfloat));
	mppf->Psi    	= (dfloat*) calloc(                               Ntotal,sizeof(dfloat));
	mppf->Rho    	= (dfloat*) calloc(                 							Ntotal,sizeof(dfloat));
	mppf->Mu     	= (dfloat*) calloc(                 							Ntotal,sizeof(dfloat));
	
	// rhs storage
	mppf->rhsU    = (dfloat*) calloc(Ntotal,sizeof(dfloat));
	mppf->rhsV    = (dfloat*) calloc(Ntotal,sizeof(dfloat));
	mppf->rhsW    = (dfloat*) calloc(Ntotal,sizeof(dfloat));
	mppf->rhsP    = (dfloat*) calloc(Ntotal,sizeof(dfloat));
	mppf->rhsPhi  = (dfloat*) calloc(Ntotal,sizeof(dfloat));

	//additional field storage
	mppf->NU      = (dfloat*) calloc(mppf->NVfields*(mppf->Nstages+1)*Ntotal,sizeof(dfloat));
	mppf->LU      = (dfloat*) calloc(mppf->NVfields*(mppf->Nstages+1)*Ntotal,sizeof(dfloat));
	mppf->GP      = (dfloat*) calloc(mppf->NVfields*(mppf->Nstages+1)*Ntotal,sizeof(dfloat));
	mppf->NPhi    = (dfloat*) calloc(               (mppf->Nstages+1)*Ntotal,sizeof(dfloat));

	mppf->GU      = (dfloat*) calloc(mppf->NVfields*Ntotal*4,sizeof(dfloat));

	mppf->rkU     = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
	mppf->rkP     = (dfloat*) calloc(               Ntotal,sizeof(dfloat));
	mppf->PI      = (dfloat*) calloc(               Ntotal,sizeof(dfloat));

	mppf->rkNU    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
	mppf->rkLU    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
	mppf->rkGP    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));

	//plotting fields
	mppf->Vort    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
	mppf->Div     = (dfloat*) calloc(               Nlocal,sizeof(dfloat));

	//  Substepping will come  here // NOT IN USE CURRENTLY
  if(mppf->elementType==HEXAHEDRA)
    mppf->cU = (dfloat *) calloc(mppf->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
  else 
    mppf->cU = mppf->U;

  mppf->Nsubsteps = 0;
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    options.getArgs("SUBCYCLING STEPS",mppf->Nsubsteps);

  if(mppf->Nsubsteps){
    mppf->Ud    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
    mppf->Ue    = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
    mppf->resU  = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));
    mppf->rhsUd = (dfloat*) calloc(mppf->NVfields*Ntotal,sizeof(dfloat));

    if(mppf->elementType==HEXAHEDRA)
      mppf->cUd = (dfloat *) calloc(mppf->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
    else 
      mppf->cUd = mppf->U;
  }



	// No gravitational acceleration
	dfloat g[3]; g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;  
  
  // Define mean-flow
	options.getArgs("UBAR", mppf->ubar);
	options.getArgs("VBAR", mppf->vbar);
	if (mppf->dim==3)
	  options.getArgs("WBAR", mppf->wbar);
	options.getArgs("PBAR", mppf->pbar);

	// Define phase properties
  options.getArgs("VISCOSITY PHASE 1", mppf->mu1);
  options.getArgs("VISCOSITY PHASE 2", mppf->mu2);

  options.getArgs("DENSITY PHASE 1", mppf->rho1);
  options.getArgs("DENSITY PHASE 2", mppf->rho2);
  
	
	//Reynolds number defined on the phase with lower material properties
  int phase_check = 1; 
  phase_check = (mppf->mu1 > mppf->mu2) ? 0: 1; 
  if(phase_check == 0){ printf("Phase 1 viscosity has to be smaller than Phase 2 \n"); exit(EXIT_FAILURE);} 

  phase_check = (mppf->rho1 > mppf->rho2) ? 0: 1; 
  if(phase_check == 0){ printf("Phase 1 density has to be smaller than Phase 2 \n"); exit(EXIT_FAILURE);} 


	mppf->Re = mppf->ubar*mppf->rho1/mppf->mu1;

	occa::properties kernelInfo;
	kernelInfo["defines"].asObject();
	kernelInfo["includes"].asArray();
	kernelInfo["header"].asArray();
	kernelInfo["flags"].asObject();

	if(mppf->dim==3)
	  meshOccaSetup3D(mesh, options, kernelInfo);
	else
	  meshOccaSetup2D(mesh, options, kernelInfo);


	occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;

  // ADD-DEFINES
  kernelInfo["defines/" "p_pbar"]= mppf->pbar;
  kernelInfo["defines/" "p_ubar"]= mppf->ubar;
  kernelInfo["defines/" "p_vbar"]= mppf->vbar;
  kernelInfo["defines/" "p_wbar"]= mppf->wbar;

  // Multiphase
  kernelInfo["defines/" "p_mu1"] = mppf->mu1;
  kernelInfo["defines/" "p_mu2"] = mppf->mu2;
  kernelInfo["defines/" "p_rho1"] = mppf->rho1;
  kernelInfo["defines/" "p_rho2"] = mppf->rho2;

  kernelInfo["defines/" "p_NTfields"]= mppf->NTfields;
  kernelInfo["defines/" "p_NVfields"]= mppf->NVfields;
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  mppf->Nstages;
  kernelInfo["defines/" "p_SUBCYCLING"]=  mppf->Nsubsteps;


  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();

  mppf->o_U 	= mesh->device.malloc(mppf->NVfields*mppf->Nstages*Ntotal*sizeof(dfloat), mppf->U);
  mppf->o_P 	= mesh->device.malloc(               mppf->Nstages*Ntotal*sizeof(dfloat), mppf->P);
  mppf->o_Phi = mesh->device.malloc(               mppf->Nstages*Ntotal*sizeof(dfloat), mppf->Phi);
  
   for (int r=0;r<mesh->size;r++) {
    if (r==mesh->rank) {
      if (mppf->dim==2) 
        mppf->setFlowFieldKernel =  mesh->device.buildKernel(DMPPF "/okl/mppfSetFlowField2D.okl", "mppfSetFlowField2D", kernelInfo);  
      // else
      //   mppf->setFlowFieldKernel =  mesh->device.buildKernel(DMPPF "/okl/mppfSetFlowField3D.okl", "mppfSetFlowField3D", kernelInfo);  
    }
    MPI_Barrier(mesh->comm);
  }





return mppf;
}