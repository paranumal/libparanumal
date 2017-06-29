#include "ellipticTri2D.h"

// create mini-meshes for each polynomial degree
void ellipticBuildAmgLevelsTri2D(solver_t *solver){

  // build some levels
  mesh2D *mesh = solver->mesh;

  // start with all polynomial levels { just for draft }
  int Nlevels = mesh->N;
  mesh2D **meshLevels = (mesh2D**) calloc(Nlevels, sizeof(mesh2D*));

  for(int level=Nlevels-1;level>=1;--level){ // hard coded for all degrees at the moment

    // hard code degree for this level
    iint levelN = level; 

    mesh2D *meshL = meshLevels[level];
    
    // copy from original mesh (important to capture geofacs, occa device, ... WILL BREAK FOR QUADS)
    memcpy(meshL, mesh, sizeof(mesh2D*));

    // reload for new degree
    meshLoadReferenceNodesTri2D(meshL, levelN);

    // set up halo exchange info for MPI (do before connect face nodes)
    meshHaloSetup(meshL);
    
    // connect face nodes (find trace indices)
    meshConnectFaceNodes2D(meshL);

    // global nodes ( do we need this ? )
    meshParallelConnectNodes(meshL);

    // ----------------------------------------------------------------------
    // specialized for matrix-free IPDG: DrT, DsT, LIFTT, MM
    // build Dr, Ds, LIFT transposes
    dfloat *DrT = (dfloat*) calloc(meshL->Np*meshL->Np, sizeof(dfloat));
    dfloat *DsT = (dfloat*) calloc(meshL->Np*meshL->Np, sizeof(dfloat));
    dfloat *LIFTT = (dfloat*) calloc(meshL->Np*meshL->Nfaces*meshL->Nfp, sizeof(dfloat));
    
    for(iint n=0;n<meshL->Np;++n){
      for(iint m=0;m<meshL->Np;++m){
	DrT[n+m*meshL->Np] = meshL->Dr[n*meshL->Np+m];
	DsT[n+m*meshL->Np] = meshL->Ds[n*meshL->Np+m];
      }
    }
    
    for(iint n=0;n<meshL->Np;++n){
      for(iint m=0;m<meshL->Nfaces*meshL->Nfp;++m){
	LIFTT[n+m*meshL->Np] = meshL->LIFT[n*meshL->Nfp*meshL->Nfaces+m];
      }
    }

    // only set up essentials on DEVICE
    meshL->o_Dr = meshL->device.malloc(meshL->Np*meshL->Np*sizeof(dfloat),  meshL->Dr);
    meshL->o_Ds = meshL->device.malloc(meshL->Np*meshL->Np*sizeof(dfloat),  meshL->Ds);
    
    meshL->o_DrT = meshL->device.malloc(meshL->Np*meshL->Np*sizeof(dfloat), DrT);
    meshL->o_DsT = meshL->device.malloc(meshL->Np*meshL->Np*sizeof(dfloat), DsT);
    
    meshL->o_LIFT  = meshL->device.malloc(meshL->Np*meshL->Nfaces*meshL->Nfp*sizeof(dfloat), meshL->LIFT);
    meshL->o_LIFTT = meshL->device.malloc(meshL->Np*meshL->Nfaces*meshL->Nfp*sizeof(dfloat), LIFTT);

    meshL->o_vmapM = meshL->device.malloc(meshL->Nelements*meshL->Nfp*meshL->Nfaces*sizeof(iint), meshL->vmapM);
    meshL->o_vmapP = meshL->device.malloc(meshL->Nelements*meshL->Nfp*meshL->Nfaces*sizeof(iint), meshL->vmapP);
    
    meshL->o_EToB =  meshL->device.malloc(meshL->Nelements*meshL->Nfaces*sizeof(iint), meshL->EToB);
    
    meshL->o_x = meshL->device.malloc(meshL->Nelements*meshL->Np*sizeof(dfloat), meshL->x);
    meshL->o_y = meshL->device.malloc(meshL->Nelements*meshL->Np*sizeof(dfloat), meshL->y);
    
    // ----------------------------------------------------------------------
    
    free(DrT); free(DsT); free(LIFTT);
  }
  
}
