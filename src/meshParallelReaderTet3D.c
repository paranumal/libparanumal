#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "mpi.h"

#include "mesh3D.h"

/* 
   purpose: read gmsh tetrahedra mesh 
*/
mesh3D* meshParallelReaderTet3D(char *fileName){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE *fp = fopen(fileName, "r");
  int n;

  mesh3D *mesh = (mesh3D*) calloc(1, sizeof(mesh3D));

  mesh->dim = 3;
  mesh->Nverts = 4; // number of vertices per element
  mesh->Nfaces = 4;
  
  // vertices on each face
  int faceVertices[4][3] = {{0,1,2},{0,1,3},{1,2,3},{2,0,3}};
  mesh->NfaceVertices = 3;
  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices*mesh->Nfaces, sizeof(int));
  memcpy(mesh->faceVertices, faceVertices[0], 12*sizeof(int));
    
  if(fp==NULL){
    printf("meshReaderTet3D: could not load file %s\n", fileName);
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

  /* find # of tets */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  int Ntets = 0, NboundaryFaces = 0;
  for(n=0;n<mesh->Nelements;++n){
    int elementType;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==4) ++Ntets; // tet code is 4
    if(elementType==2) ++NboundaryFaces;
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  int chunk = Ntets/size;
  int remainder = Ntets - chunk*size;

  int NtetsLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  int start = rank*chunk + mymin(rank, remainder); 
  int end = start + NtetsLocal-1;
  
  /* allocate space for Element node index data */

  mesh->EToV 
    = (int*) calloc(NtetsLocal*mesh->Nverts, sizeof(int));
  mesh->elementInfo
    = (int*) calloc(NtetsLocal,sizeof(int));

  /* scan through file looking for tetrahedra elements */
  int cnt=0, bcnt = 0;
  Ntets = 0;

  mesh->boundaryInfo = (int*) calloc(NboundaryFaces*4, sizeof(int));
  for(n=0;n<mesh->Nelements;++n){
    int elementType, v1, v2, v3, v4;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==2){ // boundary face
      sscanf(buf, "%*d%*d %*d%d%*d %d%d%d", 
	     mesh->boundaryInfo+bcnt*4, &v1, &v2, &v3);
      mesh->boundaryInfo[bcnt*4+1] = v1-1;
      mesh->boundaryInfo[bcnt*4+2] = v2-1;
      mesh->boundaryInfo[bcnt*4+3] = v3-1;
      ++bcnt;
    }

    if(elementType==4){  // tet code is 4
      if(start<=Ntets && Ntets<=end){
	sscanf(buf,
	       "%*d%*d%*d %d %*d"
	       intFormat intFormat intFormat intFormat,
	       mesh->elementInfo+cnt,&v1, &v2, &v3, &v4);
	/* read vertex triplet for trianngle */
	mesh->EToV[cnt*mesh->Nverts+0] = v1-1;
	mesh->EToV[cnt*mesh->Nverts+1] = v2-1;
	mesh->EToV[cnt*mesh->Nverts+2] = v3-1;
	mesh->EToV[cnt*mesh->Nverts+3] = v4-1;
	++cnt;
      }
      ++Ntets;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  mesh->NboundaryFaces = bcnt;
  
  /* record number of found tets */
  mesh->Nelements = NtetsLocal;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  for(int e=0;e<mesh->Nelements;++e){
    for(n=0;n<mesh->Nverts;++n){
      int vid = mesh->EToV[e*mesh->Nverts+n];
      mesh->EX[e*mesh->Nverts+n] = VX[vid];
      mesh->EY[e*mesh->Nverts+n] = VY[vid];
      mesh->EZ[e*mesh->Nverts+n] = VZ[vid];
    }
  }
  
  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);
  free(VZ);

  return mesh;

}
  
