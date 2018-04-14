
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "mpi.h"

#include "mesh3D.h"

/* 
   purpose: read gmsh quadrilateral mesh 
*/
mesh_t* meshParallelReaderQuad3D(char *fileName){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE *fp = fopen(fileName, "r");
  int n;

  mesh_t *mesh = (mesh_t*) calloc(1, sizeof(mesh_t));

  mesh->dim = 3;
  mesh->Nverts = 4; // number of vertices per element
  mesh->Nfaces = 4;
  mesh->NfaceVertices = 2;
     
  iint faceVertices[4][2] = {{0,1},{1,2},{2,3},{3,0}}; 
  
  mesh->faceVertices =
    (iint*) calloc(mesh->NfaceVertices*mesh->Nfaces, sizeof(iint));
  
  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices*mesh->Nfaces*sizeof(iint));
  
  if(fp==NULL){
    printf("meshReader2D: could not load file %s\n", fileName);
    exit(0);
  }

  char buf[BUFSIZ];
  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->Nnodes));

  /* allocate space for node coordinates */
  dfloat *VX = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  dfloat *VY = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  dfloat *VZ = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  
  /* load nodes */
  for(n=0;n<mesh->Nnodes;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat dfloatFormat,
	   VX+n, VY+n, VZ+n);
  }
  
  /* look for section with Element node data */
  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->Nelements));

  /* find # of quadrilaterals */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  int Nquadrilaterals = 0;

  int NboundaryFaces = 0;
  for(n=0;n<mesh->Nelements;++n){
    iint elementType;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==1) ++NboundaryFaces;
    if(elementType==3) ++Nquadrilaterals;
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  int chunk = Nquadrilaterals/size;
  int remainder = Nquadrilaterals - chunk*size;

  int NquadrilateralsLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  int start = rank*chunk + mymin(rank, remainder); 
  int end = start + NquadrilateralsLocal-1;
  
  /* allocate space for Element node index data */

  mesh->cubeFaceNumber
    = (iint *) calloc(NquadrilateralsLocal,sizeof(iint));

  mesh->cubeDistance
    = (iint *) calloc(NquadrilateralsLocal,sizeof(iint));
  
  mesh->EToV 
    = (iint*) calloc(NquadrilateralsLocal*mesh->Nverts, 
		     sizeof(iint));

  mesh->elementInfo
    = (int*) calloc(NquadrilateralsLocal,sizeof(int));

  //allocate space to store reference vertices
  //(needed to enforce rotation)

  dfloat *refNodes = (dfloat *) calloc(6*mesh->Nverts*3,sizeof(dfloat));
  
  /* scan through file looking for quadrilateral elements */
  int cnt=0, bcnt=0;
  Nquadrilaterals = 0;

  mesh->boundaryInfo = (iint*) calloc(NboundaryFaces*3, sizeof(iint));
  for(n=0;n<mesh->Nelements;++n){
    iint elementType, v1, v2, v3, v4;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);

    if(elementType==1){ // boundary face
      sscanf(buf, "%*d%*d %*d%d%*d %d%d", 
	     mesh->boundaryInfo+bcnt*3, &v1, &v2);
      mesh->boundaryInfo[bcnt*3+1] = v1-1;
      mesh->boundaryInfo[bcnt*3+2] = v2-1;
      ++bcnt;
    }
    
    if(elementType==3){  // quadrilateral
      if(start<=Nquadrilaterals && Nquadrilaterals<=end){
	
	sscanf(buf, "%*d%*d%*d %d %*d %d%d%d%d%d%d", 
	       mesh->elementInfo+cnt,mesh->cubeFaceNumber+cnt,mesh->cubeDistance+cnt,&v1, &v2, &v3, &v4);
	++cnt;
      }
      ++Nquadrilaterals;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  mesh->NboundaryFaces = bcnt;

  /* record number of found quadrilaterals */
  mesh->Nelements = NquadrilateralsLocal;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  for(int e=0;e<mesh->Nelements;++e){
    for(n=0;n<mesh->Nverts;++n){
      mesh->EX[e*mesh->Nverts+n] = VX[mesh->EToV[e*mesh->Nverts+n]];
      mesh->EY[e*mesh->Nverts+n] = VY[mesh->EToV[e*mesh->Nverts+n]];
      mesh->EZ[e*mesh->Nverts+n] = VZ[mesh->EToV[e*mesh->Nverts+n]];
#if 0
      printf("e %d v %d %g %g %g\n",
	     e, n,
	     mesh->EX[e*mesh->Nverts+n],
	     mesh->EY[e*mesh->Nverts+n],
	     mesh->EZ[e*mesh->Nverts+n]);
#endif
    }
  }

  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);
  free(VZ);

  return mesh;

}
  
