/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "mesh.hpp"

/*
   purpose: read gmsh quadrilateral mesh
*/
void meshQuad3D::ParallelReader(const char *fileName){

  FILE *fp = fopen(fileName, "r");

  char *status;

  dim = 3;
  Nverts = 4; // number of vertices per element
  Nfaces = 4;
  NfaceVertices = 2;

  int faceVertices_[4][2] = {{0,1},{1,2},{2,3},{3,0}};

  faceVertices =
    (int*) calloc(NfaceVertices*Nfaces, sizeof(int));

  memcpy(faceVertices, faceVertices_[0], NfaceVertices*Nfaces*sizeof(int));

  if(fp==NULL){
    stringstream ss;
    ss << "Cannot open file: " << fileName;
    LIBP_ABORT(ss.str())
  }

  char buf[BUFSIZ];
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &(Nnodes));

  /* allocate space for node coordinates */
  dfloat *VX = (dfloat*) calloc(Nnodes, sizeof(dfloat));
  dfloat *VY = (dfloat*) calloc(Nnodes, sizeof(dfloat));
  dfloat *VZ = (dfloat*) calloc(Nnodes, sizeof(dfloat));

  /* load nodes */
  for(int n=0;n<Nnodes;++n){
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat dfloatFormat,
	   VX+n, VY+n, VZ+n);
  }

  /* look for section with Element node data */
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(Nelements));

  /* find # of quadrilaterals */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  int Nquadrilaterals = 0;

  int NboundaryFaces = 0;
  for(int n=0;n<Nelements;++n){
    int elementType;
    status = fgets(buf, BUFSIZ, fp);
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

  EToV
    = (hlong*) calloc(NquadrilateralsLocal*Nverts,
		     sizeof(hlong));

  elementInfo
    = (hlong*) calloc(NquadrilateralsLocal,sizeof(hlong));

  /* scan through file looking for quadrilateral elements */
  int cnt=0, bcnt=0;
  Nquadrilaterals = 0;

  boundaryInfo = (hlong*) calloc(NboundaryFaces*3, sizeof(hlong));
  for(int n=0;n<Nelements;++n){
    int elementType;
    hlong v1, v2, v3, v4;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);

    if(elementType==1){ // boundary face
      sscanf(buf, "%*d%*d %*d" hlongFormat "%*d " hlongFormat hlongFormat,
	     boundaryInfo+bcnt*3, &v1, &v2);
      boundaryInfo[bcnt*3+1] = v1-1;
      boundaryInfo[bcnt*3+2] = v2-1;
      ++bcnt;
    }

    if(elementType==3){  // quadrilateral
      if(start<=Nquadrilaterals && Nquadrilaterals<=end){
        sscanf(buf, "%*d%*d%*d " hlongFormat " %*d" hlongFormat hlongFormat hlongFormat hlongFormat,
               elementInfo+cnt, &v1, &v2, &v3, &v4);

#if 0
	// check orientation
	dfloat xe1 = VX[v1-1], xe2 = VX[v2-1], xe4 = VX[v4-1];
	dfloat ye1 = VY[v1-1], ye2 = VY[v2-1], ye4 = VY[v4-1];
	dfloat J = 0.25*((xe2-xe1)*(ye4-ye1) - (xe4-xe1)*(ye2-ye1));
	if(J<0){
	  int v4tmp = v4;
	  v4 = v2;
	  v2 = v4tmp;
	  printf("unwarping element\n");
	}
#endif

	/* read vertex triplet for trianngle */
	EToV[cnt*Nverts+0] = v1-1;
	EToV[cnt*Nverts+1] = v2-1;
	EToV[cnt*Nverts+2] = v3-1;
	EToV[cnt*Nverts+3] = v4-1;
	++cnt;
      }
      ++Nquadrilaterals;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  NboundaryFaces = bcnt;

  /* record number of found quadrilaterals */
  Nelements = NquadrilateralsLocal;

  /* collect vertices for each element */
  EX = (dfloat*) calloc(Nverts*Nelements, sizeof(dfloat));
  EY = (dfloat*) calloc(Nverts*Nelements, sizeof(dfloat));
  EZ = (dfloat*) calloc(Nverts*Nelements, sizeof(dfloat));
  for(int e=0;e<Nelements;++e){
    for(int n=0;n<Nverts;++n){
      EX[e*Nverts+n] = VX[EToV[e*Nverts+n]];
      EY[e*Nverts+n] = VY[EToV[e*Nverts+n]];
      EZ[e*Nverts+n] = VZ[EToV[e*Nverts+n]];
#if 0
      printf("e %d v %d %g %g %g\n",
	     e, n,
	     EX[e*Nverts+n],
	     EY[e*Nverts+n],
	     EZ[e*Nverts+n]);
#endif
    }
  }

  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);
  free(VZ);
}

