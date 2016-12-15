
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mesh3D.h"

/* 
   purpose: read gmsh tetrahedron mesh 
*/
mesh3D* meshReader3D(char *fileName){

  FILE *fp = fopen(fileName, "r");
  int n;

  mesh3D *mesh = (mesh3D*) calloc(1, sizeof(mesh3D));

  mesh->Nverts = 4; // number of vertices per element
  mesh->Nfaces = 4;

  if(fp==NULL){
    printf("meshReader3D: could not load file %s\n", fileName);
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
    sscanf(buf, "%*d" dfloatFormat dfloatFormat dfloatFormat, VX+n, VY+n, VZ+n);
  }
  
  /* look for section with Element node data */
  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->Nelements));

  /* allocate space for Element node index data */
  mesh->EToV 
    = (iint*) calloc(mesh->Nelements*mesh->Nverts, 
		     sizeof(iint));

  /* scan through file looking for tetrahedra elements */
  int Ntets = 0;
  for(n=0;n<mesh->Nelements;++n){
    iint elementType, v1, v2, v3, v4;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==4){  // triangle
      /* not robust scanning to skip element index, element tag count, element tags */
      sscanf(buf, "%*d%*d%*d%*d%*d" iintFormat iintFormat iintFormat iintFormat,
	     &v1, &v2, &v3, &v4);
      /* read vertex triplet for trianngle */
      mesh->EToV[Ntets*mesh->Nverts+0] = v1-1;
      mesh->EToV[Ntets*mesh->Nverts+1] = v2-1;
      mesh->EToV[Ntets*mesh->Nverts+2] = v3-1;
      mesh->EToV[Ntets*mesh->Nverts+3] = v4-1;
      ++Ntets;
    }
  }
  fclose(fp);

  /* record number of found tets */
  mesh->Nelements = Ntets;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  for(int e=0;e<mesh->Nelements;++e){
    for(n=0;n<mesh->Nverts;++n){
      mesh->EX[e*mesh->Nverts+n] = VX[mesh->EToV[e*mesh->Nverts+n]];
      mesh->EY[e*mesh->Nverts+n] = VY[mesh->EToV[e*mesh->Nverts+n]];
      mesh->EZ[e*mesh->Nverts+n] = VZ[mesh->EToV[e*mesh->Nverts+n]];
    }
  }
  
  /* release VX, VY, and VZ (these are too big to keep) */
  free(VX);
  free(VY);
  free(VZ);

  return mesh;

}
  
