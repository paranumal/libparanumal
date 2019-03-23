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
   purpose: read gmsh triangle mesh
*/
void meshTri2D::ParallelReader(const char *fileName){

  FILE *fp = fopen(fileName, "r");

  char *status;

  dim = 2;
  Nverts = 3; // number of vertices per element
  Nfaces = 3;
  NfaceVertices = 2;

  /* vertices on each face */
  int faceVertices_[4][2] = {{0,1},{1,2},{2,0}};

  faceVertices =
    (int*) calloc(NfaceVertices*Nfaces, sizeof(int));

  memcpy(faceVertices, faceVertices_[0], NfaceVertices*Nfaces*sizeof(int));

  if(fp==NULL){
    stringstream ss;
    ss << "Cannot open file: " << fileName;
    LIBP_ABORT(ss.str())
  }

  char buf[BUFSIZ];


  // look for Nodes section
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &(Nnodes));

  /* allocate space for node coordinates */
  dfloat *VX = (dfloat*) calloc(Nnodes, sizeof(dfloat));
  dfloat *VY = (dfloat*) calloc(Nnodes, sizeof(dfloat));

  /* load nodes */
  for(hlong n=0;n<Nnodes;++n){
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat, VX+n, VY+n);
  }

  /* look for section with Element node data */
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of elements in mesh */
  hlong Nelements;
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &Nelements);

  /* find # of triangles */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  hlong Ntriangles = 0;
  hlong NboundaryFaces = 0;
  for(hlong n=0;n<Nelements;++n){
    int elementType;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==1) ++NboundaryFaces;
    if(elementType==2) ++Ntriangles;
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  hlong chunk = (hlong) Ntriangles/size;
  int remainder = (int) (Ntriangles - chunk*size);

  hlong NtrianglesLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  hlong start = rank*chunk + mymin(rank, remainder);
  hlong end   = start + NtrianglesLocal-1;

  /* allocate space for Element node index data */

  EToV
    = (hlong*) calloc(NtrianglesLocal*Nverts,
                     sizeof(hlong));
  elementInfo
    = (hlong*) calloc(NtrianglesLocal,sizeof(hlong));

  /* scan through file looking for triangle elements */
  hlong cnt=0, bcnt=0;
  Ntriangles = 0;

  boundaryInfo = (hlong*) calloc(NboundaryFaces*3, sizeof(hlong));
  for(hlong n=0;n<Nelements;++n){
    int elementType;
    hlong v1, v2, v3;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==1){ // boundary face
      sscanf(buf, "%*d%*d %*d" hlongFormat "%*d" hlongFormat hlongFormat,
             boundaryInfo+bcnt*3, &v1, &v2);
      boundaryInfo[bcnt*3+1] = v1-1;
      boundaryInfo[bcnt*3+2] = v2-1;
      ++bcnt;
    }
    if(elementType==2){  // triangle
      if(start<=Ntriangles && Ntriangles<=end){
        sscanf(buf, "%*d%*d%*d " hlongFormat " %*d" hlongFormat hlongFormat hlongFormat,
               elementInfo+cnt, &v1, &v2, &v3);

        // check orientation
        dfloat xe1 = VX[v1-1], xe2 = VX[v2-1], xe3 = VX[v3-1];
        dfloat ye1 = VY[v1-1], ye2 = VY[v2-1], ye3 = VY[v3-1];
        dfloat J = 0.25*((xe2-xe1)*(ye3-ye1) - (xe3-xe1)*(ye2-ye1));
        if(J<0){
          hlong v3tmp = v3;
          v3 = v2;
          v2 = v3tmp;
          //      printf("unwarping element\n");
        }

        /* read vertex triplet for trianngle */
        EToV[cnt*Nverts+0] = v1-1;
        EToV[cnt*Nverts+1] = v2-1;
        EToV[cnt*Nverts+2] = v3-1;

        ++cnt;
      }
      ++Ntriangles;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  NboundaryFaces = bcnt;

  /* record number of found triangles */
  Nelements = (dlong) NtrianglesLocal;

  /* collect vertices for each element */
  EX = (dfloat*) calloc(Nverts*Nelements, sizeof(dfloat));
  EY = (dfloat*) calloc(Nverts*Nelements, sizeof(dfloat));
  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<Nverts;++n){
      EX[e*Nverts+n] = VX[EToV[e*Nverts+n]];
      EY[e*Nverts+n] = VY[EToV[e*Nverts+n]];
    }
  }

  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);
}

