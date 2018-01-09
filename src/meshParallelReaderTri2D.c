
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "mpi.h"

#include "mesh2D.h"

/*
   purpose: read gmsh triangle mesh
*/
mesh2D* meshParallelReaderTri2D(char *fileName){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE *fp = fopen(fileName, "r");
  int n;

  char *status;

  mesh2D *mesh = (mesh2D*) calloc(1, sizeof(mesh2D));

  mesh->dim = 2;
  mesh->Nverts = 3; // number of vertices per element
  mesh->Nfaces = 3;
  mesh->NfaceVertices = 2;

  /* vertices on each face */
  iint faceVertices[4][2] = {{0,1},{1,2},{2,0}};

  mesh->faceVertices =
    (iint*) calloc(mesh->NfaceVertices*mesh->Nfaces, sizeof(iint));

  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices*mesh->Nfaces*sizeof(iint));

  if(fp==NULL){
    printf("meshParallelReaderTri2D: could not load file %s\n", fileName);
    exit(0);
  }

  char buf[BUFSIZ];


  // look for Nodes section
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->Nnodes));

  /* allocate space for node coordinates */
  dfloat *VX = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  dfloat *VY = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));

  /* load nodes */
  for(n=0;n<mesh->Nnodes;++n){
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat, VX+n, VY+n);
  }

  /* look for section with Element node data */
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->Nelements));

  /* find # of triangles */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  int Ntriangles = 0;
  int NboundaryFaces = 0;
  for(n=0;n<mesh->Nelements;++n){
    iint elementType;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==1) ++NboundaryFaces;
    if(elementType==2) ++Ntriangles;
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  int chunk = Ntriangles/size;
  int remainder = Ntriangles - chunk*size;

  int NtrianglesLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  int start = rank*chunk + mymin(rank, remainder);
  int end = start + NtrianglesLocal-1;

  /* allocate space for Element node index data */

  mesh->EToV
    = (iint*) calloc(NtrianglesLocal*mesh->Nverts,
		     sizeof(iint));
  mesh->elementInfo
    = (int*) calloc(NtrianglesLocal,sizeof(int));

  /* scan through file looking for triangle elements */
  int cnt=0, bcnt=0;
  Ntriangles = 0;

  mesh->boundaryInfo = (iint*) calloc(NboundaryFaces*3, sizeof(iint));
  for(n=0;n<mesh->Nelements;++n){
    iint elementType, v1, v2, v3;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==1){ // boundary face
      sscanf(buf, "%*d%*d %*d%d%*d %d%d",
	     mesh->boundaryInfo+bcnt*3, &v1, &v2);
      mesh->boundaryInfo[bcnt*3+1] = v1-1;
      mesh->boundaryInfo[bcnt*3+2] = v2-1;
      ++bcnt;
    }
    if(elementType==2){  // triangle
      if(start<=Ntriangles && Ntriangles<=end){
	sscanf(buf, "%*d%*d%*d %d %*d %d%d%d",
	      mesh->elementInfo+cnt, &v1, &v2, &v3);

	// check orientation
	dfloat xe1 = VX[v1-1], xe2 = VX[v2-1], xe3 = VX[v3-1];
	dfloat ye1 = VY[v1-1], ye2 = VY[v2-1], ye3 = VY[v3-1];
	dfloat J = 0.25*((xe2-xe1)*(ye3-ye1) - (xe3-xe1)*(ye2-ye1));
	if(J<0){
	  iint v3tmp = v3;
	  v3 = v2;
	  v2 = v3tmp;
	  //	  printf("unwarping element\n");
	}

	/* read vertex triplet for trianngle */
	mesh->EToV[cnt*mesh->Nverts+0] = v1-1;
	mesh->EToV[cnt*mesh->Nverts+1] = v2-1;
	mesh->EToV[cnt*mesh->Nverts+2] = v3-1;

	++cnt;
      }
      ++Ntriangles;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  mesh->NboundaryFaces = bcnt;

  /* record number of found triangles */
  mesh->Nelements = NtrianglesLocal;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  for(int e=0;e<mesh->Nelements;++e){
    for(n=0;n<mesh->Nverts;++n){
      mesh->EX[e*mesh->Nverts+n] = VX[mesh->EToV[e*mesh->Nverts+n]];
      mesh->EY[e*mesh->Nverts+n] = VY[mesh->EToV[e*mesh->Nverts+n]];
    }
  }

  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);

  return mesh;

}

