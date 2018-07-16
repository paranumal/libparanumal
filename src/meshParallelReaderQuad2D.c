
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "mpi.h"

#include "mesh2D.h"

/* 
   purpose: read gmsh quadrilateral mesh 
*/
mesh2D* meshParallelReaderQuad2D(char *fileName){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE *fp = fopen(fileName, "r");

  char *status;

  mesh2D *mesh = (mesh2D*) calloc(1, sizeof(mesh2D));

  mesh->rank = rank;
  mesh->size = size;

  
  mesh->dim = 2;
  mesh->Nverts = 4; // number of vertices per element
  mesh->Nfaces = 4;
  mesh->NfaceVertices = 2;
     
  int faceVertices[4][2] = {{0,1},{1,2},{2,3},{3,0}}; 
  
  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices*mesh->Nfaces, sizeof(int));
  
  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices*mesh->Nfaces*sizeof(int));
  
  if(fp==NULL){
    printf("meshParallelReaderQuad2D: could not load file %s\n", fileName);
    exit(0);
  }

  char buf[BUFSIZ];
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &(mesh->Nnodes));

  /* allocate space for node coordinates */
  dfloat *VX = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  dfloat *VY = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));

  /* load nodes */
  for(hlong n=0;n<mesh->Nnodes;++n){
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat, VX+n, VY+n);
  }
  
  /* look for section with Element node data */
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  hlong Nelements;
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &Nelements);

  /* find # of quadrilaterals */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  hlong Nquadrilaterals = 0;
  hlong NboundaryFaces = 0;

  for(hlong n=0;n<Nelements;++n){
    int elementType;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==1) ++NboundaryFaces;
    if(elementType==3) ++Nquadrilaterals;
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  hlong chunk = (hlong) Nquadrilaterals/size;
  int remainder = (int) (Nquadrilaterals - chunk*size);

  hlong NquadrilateralsLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  hlong start = rank*chunk + mymin(rank, remainder); 
  hlong end = start + NquadrilateralsLocal-1;
  
  /* allocate space for Element node index data */

  mesh->EToV 
    = (hlong*) calloc(NquadrilateralsLocal*mesh->Nverts, 
                     sizeof(hlong));

  mesh->elementInfo
    = (int*) calloc(NquadrilateralsLocal,sizeof(int));

  /* scan through file looking for quadrilateral elements */
  hlong cnt=0, bcnt=0;
  Nquadrilaterals = 0;

  mesh->boundaryInfo = (hlong*) calloc(NboundaryFaces*3, sizeof(hlong));
  for(hlong n=0;n<Nelements;++n){
    int elementType; 
    hlong v1, v2, v3, v4;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);

    if(elementType==1){ // boundary face
      sscanf(buf, "%*d%*d %*d" hlongFormat "%*d" hlongFormat hlongFormat,
             mesh->boundaryInfo+bcnt*3, &v1, &v2);
      mesh->boundaryInfo[bcnt*3+1] = v1-1;
      mesh->boundaryInfo[bcnt*3+2] = v2-1;
      ++bcnt;
    }
    
    if(elementType==3){  // quadrilateral
      if(start<=Nquadrilaterals && Nquadrilaterals<=end){
        sscanf(buf, "%*d%*d%*d %d %*d" hlongFormat hlongFormat hlongFormat hlongFormat, 
               mesh->elementInfo+cnt, &v1, &v2, &v3, &v4);

        // check orientation
        dfloat xe1 = VX[v1-1], xe2 = VX[v2-1], xe4 = VX[v4-1];
        dfloat ye1 = VY[v1-1], ye2 = VY[v2-1], ye4 = VY[v4-1];
        dfloat J = 0.25*((xe2-xe1)*(ye4-ye1) - (xe4-xe1)*(ye2-ye1));
        if(J<0){
          hlong v4tmp = v4;
          v4 = v2;
          v2 = v4tmp;
          //      printf("unwarping element\n");
        }
        
        /* read vertex triplet for trianngle */
        mesh->EToV[cnt*mesh->Nverts+0] = v1-1;
        mesh->EToV[cnt*mesh->Nverts+1] = v2-1;
        mesh->EToV[cnt*mesh->Nverts+2] = v3-1;
        mesh->EToV[cnt*mesh->Nverts+3] = v4-1;
        ++cnt;
      }
      ++Nquadrilaterals;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  mesh->NboundaryFaces = bcnt;
  
  /* record number of found quadrilaterals */
  mesh->Nelements = (dlong) NquadrilateralsLocal;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Nverts;++n){
      mesh->EX[e*mesh->Nverts+n] = VX[mesh->EToV[e*mesh->Nverts+n]];
      mesh->EY[e*mesh->Nverts+n] = VY[mesh->EToV[e*mesh->Nverts+n]];
    }
  }

  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);

  return mesh;

}
  
