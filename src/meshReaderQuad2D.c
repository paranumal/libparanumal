
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mesh2D.h"

/* 
   purpose: read gmsh quadrilateral mesh 
*/
mesh2D* meshReaderQuad2D(char *fileName){

  FILE *fp = fopen(fileName, "r");
  int n;

  mesh2D *mesh = (mesh2D*) calloc(1, sizeof(mesh2D));

  mesh->Nverts = 4; // number of vertices per element
  mesh->Nfaces = 4;

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

  /* load nodes */
  for(n=0;n<mesh->Nnodes;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%f%f", VX+n, VY+n);
  }
  
  /* look for section with Element node data */
  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->Nelements));

  /* allocate space for Element node index data */
  mesh->EToV = (iint*) calloc(mesh->Nelements*mesh->Nverts, sizeof(iint));

  /* scan through file looking for quadrilateral elements */
  int Nquadrilaterals = 0;
  for(n=0;n<mesh->Nelements;++n){
    iint elementType, v1, v2, v3, v4;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==3){  // quadrilateral
      /* not robust scanning to skip element index, element tag count, element tags */
      sscanf(buf, "%*d%*d%*d%*d%*d %d%d%d%d", 
	     &v1, &v2, &v3, &v4);
      /* read vertex quartet for quadrilateral */
      mesh->EToV[Nquadrilaterals*mesh->Nverts+0] = v1-1;
      mesh->EToV[Nquadrilaterals*mesh->Nverts+1] = v2-1;
      mesh->EToV[Nquadrilaterals*mesh->Nverts+2] = v3-1;
      mesh->EToV[Nquadrilaterals*mesh->Nverts+3] = v4-1;
      ++Nquadrilaterals;
    }
  }
  fclose(fp);

  /* record number of found quadrilaterals */
  mesh->Nelements = Nquadrilaterals;

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
  
