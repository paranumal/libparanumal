#include "ellipticTet3D.h"

void BuildLocalPatchAx(solver_t *solver, mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong eM, dfloat *A);

void ellipticBuildJacobiIpdgTet3D(solver_t *solver, mesh3D *mesh, int basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   int *BCType, dfloat **invDiagA,
                                   const char *options){

  if(!basis) { // default to degree N Lagrange basis
    basisNp = mesh->Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(int n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Nfp*mesh->Nfp,sizeof(dfloat));
  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Nfp;n++) {
      int fn = mesh->faceNodes[f*mesh->Nfp+n];

      for (int m=0;m<mesh->Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<mesh->Np;i++){
          MSnm += mesh->MM[fn+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
        }

        MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Nfp]  = MSnm;
      }
    }
  }

  dlong diagNnum = basisNp*mesh->Nelements;

  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  // loop over all elements
  for(dlong eM=0;eM<mesh->Nelements;++eM){
    //build the patch A matrix for this element
    BuildLocalPatchAx(solver, mesh, basis, tau, lambda, BCType, MS, eM, patchA);

    // compute the diagonal entries
    for(int n=0;n<mesh->Np;++n){
      (*invDiagA)[eM*mesh->Np + n] = 1./patchA[n*mesh->Np+n]; //store the inverse diagonal entry
    }
  }

  free(patchA);
  free(MS);
}
