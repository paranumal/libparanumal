#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "mpi.h"

#include "mesh3D.h"

/* 
   purpose: read gmsh hexrahedra mesh 
*/
mesh3D* meshParallelReaderHex3D(char *fileName){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE *fp = fopen(fileName, "r");
  int n;

  mesh3D *mesh = (mesh3D*) calloc(1, sizeof(mesh3D));

  mesh->Nverts = 8; // number of vertices per element
  mesh->Nfaces = 6;
  mesh->NfaceVertices = 4;
  
  // vertices on each face
  iint faceVertices[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},{3,4,7,6},{3,0,4,7},{4,5,6,7}};

  mesh->faceVertices =
    (iint*) calloc(mesh->NfaceVertices*mesh->Nfaces, sizeof(iint));

  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices*mesh->Nfaces*sizeof(iint));
    
  if(fp==NULL){
    printf("meshReaderHex3D: could not load file %s\n", fileName);
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

  /* find # of hexes */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  int Nhexes = 0, NboundaryFaces = 0;
  for(n=0;n<mesh->Nelements;++n){
    iint elementType;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==5) ++Nhexes; // hex code is 5
    if(elementType==3) ++NboundaryFaces; // quad codes is 3
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  int chunk = Nhexes/size;
  int remainder = Nhexes - chunk*size;

  int NhexesLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  int start = rank*chunk + mymin(rank, remainder); 
  int end = start + NhexesLocal-1;
  
  /* allocate space for Element node index data */

  mesh->EToV 
    = (iint*) calloc(NhexesLocal*mesh->Nverts, sizeof(iint));

  /* scan through file looking for hexrahedra elements */
  int cnt=0, bcnt = 0;
  Nhexes = 0;

  mesh->boundaryInfo = (iint*) calloc(NboundaryFaces*(mesh->NfaceVertices+1), sizeof(iint));
  for(n=0;n<mesh->Nelements;++n){
    iint elementType, v1, v2, v3, v4, v5, v6, v7, v8;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==3){ // quad boundary face
      iint base = (mesh->NfaceVertices+1)*bcnt;
      sscanf(buf, "%*d%*d %*d%d%*d %d%d%d%d", 
	     mesh->boundaryInfo+base, &v1, &v2, &v3, &v4);

      mesh->boundaryInfo[base+1] = v1-1;
      mesh->boundaryInfo[base+2] = v2-1;
      mesh->boundaryInfo[base+3] = v3-1;
      mesh->boundaryInfo[base+4] = v4-1;
      ++bcnt;
    }

    if(elementType==5){  // hex code is 5
      if(start<=Nhexes && Nhexes<=end){
	sscanf(buf,
	       "%*d%*d%*d%*d%*d"
	       iintFormat iintFormat iintFormat iintFormat iintFormat iintFormat iintFormat iintFormat,
	       &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8);

	mesh->EToV[cnt*mesh->Nverts+0] = v1-1;
	mesh->EToV[cnt*mesh->Nverts+1] = v2-1;
	mesh->EToV[cnt*mesh->Nverts+2] = v3-1;
	mesh->EToV[cnt*mesh->Nverts+3] = v4-1;
	mesh->EToV[cnt*mesh->Nverts+4] = v5-1;
	mesh->EToV[cnt*mesh->Nverts+5] = v6-1;
	mesh->EToV[cnt*mesh->Nverts+6] = v7-1;
	mesh->EToV[cnt*mesh->Nverts+7] = v8-1;

	//	printf("%d: %d,%d,%d,%d %d,%d,%d,%d", cnt, v1-1, v2-1,v3-1,v4-1,v5-1,v6-1,v7-1,v8-1);
	
	++cnt;
      }
      ++Nhexes;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  mesh->NboundaryFaces = bcnt;
  
  /* record number of found hexes */
  mesh->Nelements = NhexesLocal;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  for(int e=0;e<mesh->Nelements;++e){
    for(n=0;n<mesh->Nverts;++n){
      iint vid = mesh->EToV[e*mesh->Nverts+n];
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
  
