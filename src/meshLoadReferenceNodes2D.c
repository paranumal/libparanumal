
#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshLoadReferenceNodes2D(mesh2D *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, "nodes/triangleN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  char buf[BUFSIZ];
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Ncheck;
  sscanf(buf, "%d", &Ncheck);
  if(Ncheck != N) printf("bu55er - wrong data file\n");
  mesh->N = N;
  mesh->Nfp = N+1;

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Npcheck;
  sscanf(buf, "%d", &Npcheck);
  mesh->Np = Npcheck;

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->r = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  mesh->s = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat, mesh->r+n, mesh->s+n);
  }

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->Dr = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->Dr+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->Ds = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->Ds+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment

  fgets(buf, BUFSIZ, fp); // read comment

  mesh->faceNodes = (iint*) calloc(mesh->Nfp*mesh->Nfaces, sizeof(iint));
  for(int n=0;n<mesh->Nfaces*mesh->Nfp;++n){
    fscanf(fp, iintFormat, mesh->faceNodes+n);
  }
  fgets(buf, BUFSIZ, fp);

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->LIFT = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Nfaces*mesh->Nfp*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->LIFT+n);
  }
  fgets(buf, BUFSIZ, fp);

  // read number of plot nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, iintFormat, &(mesh->plotNp));

  // read plot node coordinates (hard code triangles)
  mesh->plotR = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  mesh->plotS = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNp;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat, mesh->plotR+n, mesh->plotS+n);
  }
  
  // read plot interpolation matrix
  mesh->plotInterp = (dfloat*) calloc(mesh->plotNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNp;++n){
    for(int m=0;m<mesh->Np;++m){
      fscanf(fp, dfloatFormat, mesh->plotInterp+n*mesh->Np+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read number of elements in plot node triangulation
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, iintFormat, &(mesh->plotNelements));

  // read number of vertices per plot element
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, iintFormat, &(mesh->plotNverts));
  
  // build and read in plot node triangulation
  mesh->plotEToV = (iint*) calloc(mesh->plotNelements*mesh->plotNverts, sizeof(iint));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNelements;++n){
    for(int m=0;m<mesh->plotNverts;++m){
      fscanf(fp, iintFormat, mesh->plotEToV+m + mesh->plotNverts*n);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read number of volume cubature nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, iintFormat, &(mesh->cubNp));

  // read cub nodes and weights
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->cubr = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  mesh->cubs = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  for(int n=0;n<mesh->cubNp;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat, mesh->cubr+n, mesh->cubs+n);
  } 

  // read volume cubature interpolation matrix
  mesh->cubInterp = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment

  for(int n=0;n<mesh->cubNp;++n){
    for(int m=0;m<mesh->Np;++m){
      fscanf(fp, dfloatFormat, mesh->cubInterp+n*mesh->Np+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read cubature weak 'r' differentiation matrix
  mesh->cubDrW = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubDrW+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }
  // read cubature weak 's' differentiation matrix
  mesh->cubDsW = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubDsW+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

    // read cubature projection matrix
  mesh->cubProject = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubProject+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }


  // read number of surface integration nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, iintFormat, &(mesh->intNfp));

  // read surface intergration node interpolation matrix
  mesh->intInterp 
    = (dfloat*) calloc(mesh->intNfp*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment

  for(int n=0;n<mesh->intNfp*mesh->Nfaces;++n){
    for(int m=0;m<mesh->Nfp;++m){
      fscanf(fp, dfloatFormat, mesh->intInterp+n*mesh->Nfp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read lift matrix from surface integration to volume nodes
  mesh->intLIFT = (dfloat*) calloc(mesh->intNfp*mesh->Nfaces*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->intNfp*mesh->Nfaces;++m){
      fscanf(fp, dfloatFormat, mesh->intLIFT+n*mesh->intNfp*mesh->Nfaces+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // ======================= BB data [added by JC] ===========================
#if 1
  // BB Vandermonde matrices
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->VB = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->VB+n);
  }

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->invVB = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->invVB+n);
  }

  // sparse barycentric differentiation matrices
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->D1ids = (iint*) calloc(mesh->Np*3, sizeof(iint));
  for (int n=0;n<mesh->Np*3;++n){    
    fscanf(fp, iintFormat, mesh->D1ids+n);
  }

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->D2ids = (iint*) calloc(mesh->Np*3, sizeof(iint));
  for (int n=0;n<mesh->Np*3;++n){
    fscanf(fp, iintFormat, mesh->D2ids+n);
  }
  
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->D3ids = (iint*) calloc(mesh->Np*3, sizeof(iint));
  for (int n=0;n<mesh->Np*3;++n){
    fscanf(fp, iintFormat, mesh->D3ids+n);
  }

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->Dvals = (dfloat*) calloc(mesh->Np*3, sizeof(dfloat));
  for (int n=0;n<mesh->Np*3;++n){
    fscanf(fp, dfloatFormat, mesh->Dvals+n);
  }
  
  // BB cubature projection matrices
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->VBq = (dfloat*) calloc(mesh->Np*mesh->cubNp, sizeof(dfloat));
  for (int n=0;n<mesh->Np*mesh->cubNp;++n){
    fscanf(fp, dfloatFormat, mesh->VBq+n);
  }
  
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->PBq = (dfloat*) calloc(mesh->Np*mesh->cubNp, sizeof(dfloat));
  for (int n=0;n<mesh->Np*mesh->cubNp;++n){
    fscanf(fp, dfloatFormat, mesh->PBq+n);
  }

  // BB L0 matrix values
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->L0vals = (dfloat*) calloc(mesh->Nfp*3, sizeof(dfloat));
  for (int n=0;n<mesh->Nfp*3;++n){
    fscanf(fp, dfloatFormat, mesh->L0vals+n);
  }

  // BB lift reduction (EL) matrix (sparse format)
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->max_EL_nnz = mesh->Nfp + 2;
  mesh->ELids = (iint*) calloc(mesh->Np*mesh->max_EL_nnz, sizeof(iint));
  for (int n=0;n<mesh->Np*mesh->max_EL_nnz;++n){
    fscanf(fp, iintFormat, mesh->ELids+n);
  }
  
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->ELvals = (dfloat*) calloc(mesh->Np*mesh->max_EL_nnz, sizeof(dfloat));
  for (int n=0;n<mesh->Np*mesh->max_EL_nnz;++n){
    fscanf(fp, dfloatFormat, mesh->ELvals+n);
  }
  
  // ============ end BB stuff ==================
#endif
  
#if 0
  printf("Dr: \n");
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      printf("%g ", mesh->Dr[n*mesh->Np+m]);
    }
    printf("\n");
  }

  printf("Ds: \n");
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      printf("%g ", mesh->Ds[n*mesh->Np+m]);
    }
    printf("\n");
  }
#endif
  

  fclose(fp);
}


