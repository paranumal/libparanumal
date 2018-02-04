
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
	sscanf(buf, "%*d%*d%*d %d %*d %d%d%d%d%d", 
	       mesh->elementInfo+cnt,mesh->cubeFaceNumber+cnt,&v1, &v2, &v3, &v4);
	
	// check orientation using a*(bxc) > 0
	dfloat xe1 = VX[v1-1], xe2 = VX[v2-1], xe3 = VX[v3-1];
	dfloat ye1 = VY[v1-1], ye2 = VY[v2-1], ye3 = VY[v3-1];
	dfloat ze1 = VZ[v1-1], ze2 = VZ[v2-1], ze3 = VZ[v3-1];
	dfloat J = xe1*ye2*ze3 - xe1*ze2*ye3 + ye1*ze2*xe3 - ye1*xe2*ze3 + ze1*xe2*ye3 - ze1*ye2*xe3;
	if (J<0) {
	  iint v4tmp = v4;
	  v4 = v2;
	  v2 = v4tmp;
	  printf("unwarping element\n");
	}
	
	int faceId = mesh->cubeFaceNumber[cnt];
	
	//check if we have a new reference node
	if (refNodes[faceId*12] == 0 && refNodes[faceId*12+1] == 0 && refNodes[faceId*12+2] == 0) {
	  //printf("adding ref node %d\n",faceId);
	  
	  refNodes[faceId*12] = VX[v1 - 1];
	  refNodes[faceId*12+1] = VY[v1 - 1];
	  refNodes[faceId*12+2] = VZ[v1 - 1];

	  refNodes[faceId*12+3] = VX[v2 - 1];
	  refNodes[faceId*12+4] = VY[v2 - 1];
	  refNodes[faceId*12+5] = VZ[v2 - 1];

	  refNodes[faceId*12+6] = VX[v3 - 1];
	  refNodes[faceId*12+7] = VY[v3 - 1];
	  refNodes[faceId*12+8] = VZ[v3 - 1];

	  refNodes[faceId*12+9] = VX[v4 - 1];
	  refNodes[faceId*12+10] = VY[v4 - 1];
	  refNodes[faceId*12+11] = VZ[v4 - 1];

	  /* read vertex triplet for trianngle */
	  mesh->EToV[cnt*mesh->Nverts+0] = v1-1;
	  mesh->EToV[cnt*mesh->Nverts+1] = v2-1;
	  mesh->EToV[cnt*mesh->Nverts+2] = v3-1;
	  mesh->EToV[cnt*mesh->Nverts+3] = v4-1;
	  ++cnt;
	}
	//otherwise, check orientation against reference
	//maximize ((a1 + a2)x(b1+a2))*((a1+a2)x(a1+b2))
	else {
	  dfloat xe1 = VX[v1-1], xe2 = VX[v2-1], xe3 = VX[v3-1], xe4 = VX[v4-1];
	  dfloat ye1 = VY[v1-1], ye2 = VY[v2-1], ye3 = VY[v3-1], ye4 = VY[v4-1];
	  dfloat ze1 = VZ[v1-1], ze2 = VZ[v2-1], ze3 = VZ[v3-1], ze4 = VZ[v4-1];

	  dfloat xr1 = refNodes[faceId*12], xr2 = refNodes[faceId*12+3],
	         xr3 = refNodes[faceId*12+6], xr4 = refNodes[faceId*12+9];
	  dfloat yr1 = refNodes[faceId*12+1], yr2 = refNodes[faceId*12+4],
	         yr3 = refNodes[faceId*12+7], yr4 = refNodes[faceId*12+10];
	  dfloat zr1 = refNodes[faceId*12+2], zr2 = refNodes[faceId*12+5],
         	 zr3 = refNodes[faceId*12+8], zr4 = refNodes[faceId*12+11];
	  int maxNode = -1;
	  dfloat maxVal = -1;

	  /*printf("refx: %lf %lf %lf %lf \n",xr1,xr2,xr3,xr4);
	  printf("refy: %lf %lf %lf %lf \n",yr1,yr2,yr3,yr4);
	  printf("refz: %lf %lf %lf %lf \n",zr1,zr2,zr3,zr4);

	  printf("elemx: %lf %lf %lf %lf \n",xe1,xe2,xe3,xe4);
	  printf("elemy: %lf %lf %lf %lf \n",ye1,ye2,ye3,ye4);
	  printf("elemz: %lf %lf %lf %lf \n",ze1,ze2,ze3,ze4);*/
	  
	  //test first edge pairing
	  
	  //find elements of cross product
	  dfloat cross1x = (ye1+yr1)*(ze2+zr1) - (ze1+zr1)*(ye2+yr1);
	  dfloat cross1y = (ze1+zr1)*(xr1+xe2) - (xe1+xr1)*(ze2+zr1);
	  dfloat cross1z = (xe1+xr1)*(ye2+yr1) - (ye1+yr1)*(xe2+xr1);
	  
	  dfloat cross2x = (ye1+yr1)*(zr2+ze1) - (zr1+ze1)*(yr2+ye1);
	  dfloat cross2y = (ze1+zr1)*(xr2+xe1) - (xr1+xe1)*(zr2+ze1);
	  dfloat cross2z = (xe1+xr1)*(yr2+ye1) - (yr1+ye1)*(xr2+xe1);

	  dfloat dotp = cross1x*cross2x + cross1y*cross2y + cross1z*cross2z;

	  //printf("%lf   ",dotp);
	  
	  maxVal = mymax(maxVal,dotp);
	  if (maxVal == dotp) maxNode = 0;

	  //test second edge pairing

	  //find elements of cross product
	  cross1x = (ye2+yr1)*(ze3+zr1) - (ze2+zr1)*(ye3+yr1);
	  cross1y = (ze2+zr1)*(xr1+xe3) - (xe2+xr1)*(ze3+zr1);
	  cross1z = (xe2+xr1)*(ye3+yr1) - (ye2+yr1)*(xe3+xr1);
	  
	  cross2x = (ye2+yr1)*(zr2+ze2) - (zr1+ze2)*(yr2+ye2);
	  cross2y = (ze2+zr1)*(xr2+xe2) - (xr1+xe2)*(zr2+ze2);
	  cross2z = (xe2+xr1)*(yr2+ye2) - (yr1+ye2)*(xr2+xe2);

	  dotp = cross1x*cross2x + cross1y*cross2y + cross1z*cross2z;

	  //printf("%lf   ",dotp);
	  
	  maxVal = mymax(maxVal,dotp);
	  if (maxVal == dotp) maxNode = 1;

	  //test third edge pairing

	  //find elements of cross product
	  cross1x = (ye3+yr1)*(ze4+zr1) - (ze3+zr1)*(ye4+yr1);
	  cross1y = (ze3+zr1)*(xr1+xe4) - (xe3+xr1)*(ze4+zr1);
	  cross1z = (xe3+xr1)*(ye4+yr1) - (ye3+yr1)*(xe4+xr1);
	  
	  cross2x = (ye3+yr1)*(zr2+ze3) - (zr1+ze3)*(yr2+ye3);
	  cross2y = (ze3+zr1)*(xr2+xe3) - (xr1+xe3)*(zr2+ze3);
	  cross2z = (xe3+xr1)*(yr2+ye3) - (yr1+ye3)*(xr2+xe3);

	  dotp = cross1x*cross2x + cross1y*cross2y + cross1z*cross2z;

	  //printf("%lf   ",dotp);
	  
	  maxVal = mymax(maxVal,dotp);
	  if (maxVal == dotp) maxNode = 2;	  

	  //test final edge pairing

	  //find elements of cross product
	  cross1x = (ye4+yr1)*(ze1+zr1) - (ze4+zr1)*(ye1+yr1);
	  cross1y = (ze4+zr1)*(xr1+xe1) - (xe4+xr1)*(ze1+zr1);
	  cross1z = (xe4+xr1)*(ye1+yr1) - (ye4+yr1)*(xe1+xr1);
	  
	  cross2x = (ye4+yr1)*(zr2+ze4) - (zr1+ze4)*(yr2+ye4);
	  cross2y = (ze4+zr1)*(xr2+xe4) - (xr1+xe4)*(zr2+ze4);
	  cross2z = (xe4+xr1)*(yr2+ye4) - (yr1+ye4)*(xr2+xe4);

	  dotp = cross1x*cross2x + cross1y*cross2y + cross1z*cross2z;

	  //printf("%lf   ",dotp);

	  maxVal = mymax(maxVal,dotp);
	  if (maxVal == dotp) maxNode = 3;

	  if (maxNode == -1) {printf("Bad element alignment. Check face numbering\n");}
	  else {
	    /* read shifted vertex triplet for trianngle */
	    mesh->EToV[cnt*mesh->Nverts+((0+maxNode)%4)] = v1-1;
	    mesh->EToV[cnt*mesh->Nverts+((1+maxNode)%4)] = v2-1;
	    mesh->EToV[cnt*mesh->Nverts+((2+maxNode)%4)] = v3-1;
	    mesh->EToV[cnt*mesh->Nverts+((3+maxNode)%4)] = v4-1;
	    ++cnt;
	  }
	}
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
  
