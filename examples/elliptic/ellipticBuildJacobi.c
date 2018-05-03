#include "elliptic.h"

void BuildLocalIpdgPatchAxTri2D(elliptic_t* elliptic, mesh_t *mesh, int basisNp, dfloat *basis, dfloat lambda, 
                        dfloat *MS, dlong eM, dfloat *A);

void BuildLocalContinuousPatchAxTri2D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda,
                                  dlong eM, dfloat *Srr, dfloat *Srs, 
                                  dfloat *Sss, dfloat *MM, dfloat *A);


void ellipticBuildJacobiTri2D(elliptic_t* elliptic, int basisNp, dfloat *basis, 
                                dfloat lambda, dfloat **invDiagA);


void ellipticBuildJacobi(elliptic_t* elliptic, int basisNp, dfloat *basis, 
                                dfloat lambda, dfloat **invDiagA) {

  switch(elliptic->elementType){
  case TRIANGLES:
    ellipticBuildJacobiTri2D(elliptic, basisNp, basis, lambda, invDiagA); break;
  case QUADRILATERALS:
    //ellipticBuildIpdgQuad2D(elliptic, basisNp, basis,lambda, A, nnzA, globalStarts); 
  break;
  case TETRAHEDRA:
    //ellipticBuildIpdgTet3D(elliptic, basisNp, basis,lambda, A, nnzA, globalStarts); 
  break;
  case HEXAHEDRA:
    //ellipticBuildIpdgHex3D(elliptic, basisNp, basis,lambda, A, nnzA, globalStarts); 
  break;
  }

}


void ellipticBuildJacobiTri2D(elliptic_t* elliptic, int basisNp, dfloat *basis, 
                                dfloat lambda, dfloat **invDiagA){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

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

  dfloat *Srr, *Srs, *Sss, *MM;
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    MM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Np;m++) {
        Srr[m+n*mesh->Np] = mesh->Srr[m+n*mesh->Np];
        Srs[m+n*mesh->Np] = mesh->Srs[m+n*mesh->Np] + mesh->Ssr[m+n*mesh->Np];
        Sss[m+n*mesh->Np] = mesh->Sss[m+n*mesh->Np];
        MM[m+n*mesh->Np] = mesh->MM[m+n*mesh->Np];
      }
    }
  }

  dlong diagNnum = basisNp*mesh->Nelements;

  dfloat *diagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));

  if(rank==0) printf("Building diagonal...");fflush(stdout);

  // loop over all elements
  #pragma omp parallel
  {
    //temp patch storage
    dfloat *patchA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

    #pragma omp for 
    for(dlong eM=0;eM<mesh->Nelements;++eM){
      //build the patch A matrix for this element
      if (options.compareArgs("DISCRETIZATION","IPDG")) {
        BuildLocalIpdgPatchAxTri2D(elliptic, mesh, basisNp, basis, lambda, MS, eM, patchA);
      } else if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
        BuildLocalContinuousPatchAxTri2D(elliptic, mesh, lambda, eM, Srr, Srs, Sss, MM, patchA);
      }

      for(int n=0;n<mesh->Np;++n) {
        diagA[eM*mesh->Np + n] = patchA[n*mesh->Np+n]; //store the diagonal entry
      }
    }
    free(patchA);
  }

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) 
    gsParallelGatherScatter(mesh->hostGsh, diagA, dfloatString, "add"); 
    
  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    (*invDiagA)[n] = 1/diagA[n];
  }

  if(rank==0) printf("done.\n");

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    free(Srr);
    free(Srs);
    free(Sss);
    free(MM);
  }

  free(diagA);
  free(MS);
}
