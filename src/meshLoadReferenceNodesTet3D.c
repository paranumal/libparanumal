
#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void meshLoadReferenceNodesTet3D(mesh3D *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/nodes/tetN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  char buf[BUFSIZ];
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Ncheck;
  sscanf(buf, "%d", &Ncheck);
  if(Ncheck != N) printf("bu55er - wrong data file\n");

  mesh->N = N;
  mesh->Np = ((N+1)*(N+2)*(N+3))/6;
  mesh->Nfp = ((N+1)*(N+2))/2;

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Npcheck;
  sscanf(buf, "%d", &Npcheck);
  mesh->Np = Npcheck;

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->r = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  mesh->s = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  mesh->t = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat dfloatFormat,
	   mesh->r+n, mesh->s+n, mesh->t+n);
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
  mesh->Dt = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->Dt+n);
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
  fgets(buf, BUFSIZ, fp); // read comment
  
  // read number of plot nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, iintFormat, &(mesh->plotNp));

  // read plot node coordinates (hard code triangles)
  mesh->plotR = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  mesh->plotS = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  mesh->plotT = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNp;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat dfloatFormat, mesh->plotR+n, mesh->plotS+n, mesh->plotT+n);
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
  printf("%s", buf);
  fgets(buf, BUFSIZ, fp);
  printf("%s", buf);
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

  fgets(buf, BUFSIZ, fp); // read comment
  printf("%s", buf);
  fgets(buf, BUFSIZ, fp); // read number of cubature points
  printf("%s", buf);
  sscanf(buf, iintFormat, &(mesh->cubNp));
  printf("cubNp = %d\n",mesh->cubNp);
    
  mesh->cubr = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  mesh->cubs = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  mesh->cubt = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  printf("%s", buf);
  for(int n=0;n<mesh->cubNp;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat dfloatFormat, mesh->cubr+n, mesh->cubs+n, mesh->cubt+n);
  }
  
  mesh->cubInterp = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->cubNp*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->cubInterp+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment  

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


   // read cubature weak 's' differentiation matrix
  mesh->cubDtW = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubDtW+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }
  fgets(buf, BUFSIZ, fp); // read comment  

  
  mesh->cubProject = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->cubNp*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->cubProject+n);
  }

  //-------------Berstein Bezier DG stuff added by NC--------------------//
  mesh->VB = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->VB+n);
  }
  fgets(buf, BUFSIZ, fp); 
  
  mesh->invVB = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->invVB+n);
  }
  fgets(buf, BUFSIZ, fp);   

  mesh->D0ids = (iint*) calloc(mesh->Np*4, sizeof(iint));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np*4;++n){
    fscanf(fp, iintFormat, mesh->D0ids+n);
  }
  fgets(buf, BUFSIZ, fp);   

  mesh->D1ids = (iint*) calloc(mesh->Np*4, sizeof(iint));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np*4;++n){
    fscanf(fp, iintFormat, mesh->D1ids+n);
  }
  fgets(buf, BUFSIZ, fp);   

  mesh->D2ids = (iint*) calloc(mesh->Np*4, sizeof(iint));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np*4;++n){
    fscanf(fp, iintFormat, mesh->D2ids+n);
  }
  fgets(buf, BUFSIZ, fp);   

  mesh->D3ids = (iint*) calloc(mesh->Np*4, sizeof(iint));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np*4;++n){
    fscanf(fp, iintFormat, mesh->D3ids+n);
  }
  fgets(buf, BUFSIZ, fp);   

  mesh->Dvals = (dfloat*) calloc(mesh->Np*4, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np*4;++n){
    fscanf(fp, dfloatFormat, mesh->Dvals+n);
  }
  fgets(buf, BUFSIZ, fp);   

  mesh->L0ids = (iint*) calloc(mesh->Nfp*7, sizeof(iint));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Nfp*7;++n){
    fscanf(fp, iintFormat, mesh->L0ids+n);
  }
  fgets(buf, BUFSIZ, fp);   

  mesh->L0vals = (dfloat*) calloc(mesh->Nfp*7, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Nfp*7;++n){
    fscanf(fp, dfloatFormat, mesh->L0vals+n);
  }
  fgets(buf, BUFSIZ, fp);   

  mesh->max_EL_nnz = mesh->Nfp+3;
  mesh->ELids = (iint*) calloc(mesh->Np*mesh->max_EL_nnz, sizeof(iint));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np*mesh->max_EL_nnz;++n){
    fscanf(fp, iintFormat, mesh->ELids+n);
  }
  fgets(buf, BUFSIZ, fp);   

  mesh->ELvals = (dfloat*) calloc(mesh->Np*mesh->max_EL_nnz, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np*mesh->max_EL_nnz;++n){
    fscanf(fp, dfloatFormat, mesh->ELvals+n);
  }
  fgets(buf, BUFSIZ, fp); 
  
#if 0
  for(int n=0;n<mesh->cubNp;++n){
    printf("rq,sq,tq = %f, %f, %f\n",mesh->cubr[n],mesh->cubs[n],mesh->cubt[n]);
  }

  printf("Vq: \n");
  for(int n=0;n<mesh->cubNp;++n){
    for(int m=0;m<mesh->Np;++m){
      printf("%f ",mesh->cubInterp[n*mesh->Np+m]);
    }
    printf("\n");
  }

  printf("Pq: \n");
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      printf("%f ",mesh->cubProject[n*mesh->cubNp+m]);
    }
    printf("\n");
  }
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

  printf("Dt: \n");
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      printf("%g ", mesh->Dt[n*mesh->Np+m]);
    }
    printf("\n");
  }
  
  printf("LIFT: \n");
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Nfaces*mesh->Nfp;++m){
      printf("%g ", mesh->LIFT[n*mesh->Nfaces*mesh->Nfp+m]);
    }
    printf("\n");
  }

  printf("faceNodes: \n");
  for(int f=0;f<mesh->Nfaces;++f){
    for(int n=0;n<mesh->Nfp;++n){
      printf("%d ", mesh->faceNodes[n+f*mesh->Nfp]);
    }
    printf("\n");
  }
#endif

  fclose(fp);
}
