#include "advectionQuad3D.h"

void advectionRunLSERKQuad3D(solver_t *solver){

  mesh_t *mesh = solver->mesh;


  occa::memory o_qOrig = mesh->device.malloc(mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat));
  
  // some sanity checks
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat J = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + JID*mesh->Np + n];
      if(J<1e-10)
	printf("J = %\g", J);
    }

#if 0
    
    int isInternal = 1;
    int cubeFace = mesh->cubeFaceNumber[e];

    for(int f=0;f<mesh->Nfaces;++f){
      int cubeFaceNeighbor = mesh->cubeFaceNumber[mesh->EToE[mesh->Nfaces*e+f]];
      if(cubeFace!=cubeFaceNeighbor)
	isInternal = 0;
    }

    if(isInternal){
      printf("EToF[%d,*] = [", e);
      for(int f=0;f<mesh->Nfaces;++f){
	printf("%d ", mesh->EToF[mesh->Nfaces*e+f]);
      }
      printf("]\n");
    }
#endif

    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->Nq;++n){
	int id = e*mesh->Nfaces*mesh->Nq + f*mesh->Nq + n;
	int vidM = mesh->vmapM[id];
	int vidP = mesh->vmapP[id];

	const dfloat dx = mesh->x[vidP]-mesh->x[vidM];
	const dfloat dy = mesh->y[vidP]-mesh->y[vidM];
	const dfloat dz = mesh->z[vidP]-mesh->z[vidM];

	if(dx*dx+dy*dy+dz*dz>1e-10)
	  printf("%lg,%lg,%lg \n",dx,dy,dz);
      }
    }
    
  }


  
  // MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  dfloat * test_q = (dfloat *) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
    
  //kernel arguments
  dfloat alpha = 0; 1./mesh->N;

  printf("alpha=%g\n", alpha);

  mesh->o_qpre.copyTo(mesh->q);
  advectionErrorNormQuad3D(mesh, 0, "icPreFilter", 0);
  
  mesh->filterKernelq0H(mesh->Nelements,
			alpha,
			mesh->o_dualProjMatrix,
			mesh->o_cubeFaceNumber,
			mesh->o_EToE,
			mesh->o_x,
			mesh->o_y,
			mesh->o_z,
			mesh->o_qpre,
			o_qOrig,
			mesh->o_qPreFilter);

  mesh->o_qPreFilter.copyTo(mesh->q);
  advectionErrorNormQuad3D(mesh, 0, "icMiddleFilter", 0);
  
  mesh->filterKernelq0V(mesh->Nelements,
			alpha,
			mesh->o_dualProjMatrix,
			mesh->o_cubeFaceNumber,
			mesh->o_EToE,
			mesh->o_x,
			mesh->o_y,
			mesh->o_z,
			mesh->o_qPreFilter,
			o_qOrig,
			mesh->o_qpre);

  mesh->o_qpre.copyTo(mesh->q);
  advectionErrorNormQuad3D(mesh, 0, "icPostFilter", 0);
  
  for(iint tstep=0;tstep < mesh->NtimeSteps;++tstep){
    for (iint Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {

      iint lev;
      for (lev=0;lev<mesh->MRABNlevels;lev++)
	if (Ntick % (1<<lev) !=0) break; //find the max lev to add to rhsq

      iint levS;
      for (levS=0;levS<mesh->MRABNlevels;levS++)
	if ((Ntick+1) % (1<<levS) !=0) break; //find the max lev to add to rhsq

      // TW: I am not sure why this pre-filter fixes the stability issue
      
      mesh->filterKernelq0H(mesh->Nelements,
			    alpha,
			    mesh->o_dualProjMatrix,
			    mesh->o_cubeFaceNumber,
			    mesh->o_EToE,
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    mesh->o_qpre,
			    o_qOrig,
			    mesh->o_qPreFilter);

      mesh->filterKernelq0V(mesh->Nelements,
			    alpha,
			    mesh->o_dualProjMatrix,
			    mesh->o_cubeFaceNumber,
			    mesh->o_EToE,
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    mesh->o_qPreFilter,
			    o_qOrig,
			    mesh->o_qpre);
      
      for (iint rk = 0; rk < mesh->Nrk; ++rk) {
	
	//synthesize actual stage time
	dfloat t = tstep*pow(2,mesh->MRABNlevels-1) + Ntick;

	// compute volume contribution to DG advection RHS
	mesh->volumePreKernel(mesh->Nelements,
			      mesh->o_vgeo,
			      mesh->o_D,
			      mesh->o_x,
			      mesh->o_y,
			      mesh->o_z,
			      mesh->o_qpre,
			      mesh->o_prerhsq);

	mesh->surfacePreKernel(mesh->Nelements,
			       mesh->o_sgeo,
			       mesh->o_LIFTT,
			       mesh->o_vmapM,
			       mesh->o_vmapP,
			       t,
			       mesh->o_x,
			       mesh->o_y,
			       mesh->o_z,
			       mesh->o_qpre,
			       mesh->o_prerhsq);

	mesh->filterKernelq0H(mesh->Nelements,
			      alpha,
			      mesh->o_dualProjMatrix,
			      mesh->o_cubeFaceNumber,
			      mesh->o_EToE,
			      mesh->o_x,
			      mesh->o_y,
			      mesh->o_z,
			      mesh->o_prerhsq,
			      o_qOrig,
			      mesh->o_qPreFilter);

	mesh->filterKernelq0V(mesh->Nelements,
			      alpha,
			      mesh->o_dualProjMatrix,
			      mesh->o_cubeFaceNumber,
			      mesh->o_EToE,
			      mesh->o_x,
			      mesh->o_y,
			      mesh->o_z,
			      mesh->o_qPreFilter,
			      o_qOrig,
			      mesh->o_prerhsq);

	
	mesh->updatePreKernel(mesh->Nelements,
			      mesh->dt,
			      mesh->rka[rk],
			      mesh->rkb[rk],
			      mesh->o_prerhsq,
			      mesh->o_resq,
			      mesh->o_qpre);
      }
    }
    
    if (tstep==0 || (tstep+1) % mesh->errorStep == 0) {
      mesh->o_qpre.copyTo(mesh->q);
      dfloat t = mesh->dt*(tstep+1)*pow(2,mesh->MRABNlevels-1);
      advectionErrorNormQuad3D(mesh, t, "foo", (tstep+1)/mesh->errorStep);
    }

  }

  
  free(recvBuffer);
  free(sendBuffer);
}
