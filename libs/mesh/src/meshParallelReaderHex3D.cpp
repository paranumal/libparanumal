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
   purpose: read gmsh hexrahedra mesh
*/
void meshHex3D::ParallelReader(const char *fileName){

  FILE *fp = fopen(fileName, "r");

  char *status;

  dim = 3;
  Nverts = 8; // number of vertices per element
  Nfaces = 6;
  NfaceVertices = 4;

  // vertices on each face
  int _faceVertices[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7}};

  faceVertices =
    (int*) calloc(NfaceVertices*Nfaces, sizeof(int));

  memcpy(faceVertices, _faceVertices[0], NfaceVertices*Nfaces*sizeof(int));

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
  for(hlong n=0;n<Nnodes;++n){
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat dfloatFormat,
           VX+n, VY+n, VZ+n);

  }

  /* look for section with Element node data */
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  hlong Nelements;
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &Nelements);

  /* find # of hexes */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  hlong Nhexes = 0, NboundaryFaces = 0;

  for(hlong n=0;n<Nelements;++n){
    int elementType;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==5) ++Nhexes; // hex code is 5
    if(elementType==3) ++NboundaryFaces; // quad codes is 3
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  hlong chunk = (hlong) Nhexes/size;
  int remainder = (int) (Nhexes - chunk*size);

  hlong NhexesLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  hlong start = rank*chunk + mymin(rank, remainder);
  hlong end = start + NhexesLocal-1;

  /* allocate space for Element node index data */

  EToV
    = (hlong*) calloc(NhexesLocal*Nverts, sizeof(hlong));

  elementInfo
    = (hlong*) calloc(NhexesLocal,sizeof(hlong));

  /* scan through file looking for hexrahedra elements */
  hlong cnt=0, bcnt=0;
  Nhexes = 0;

  boundaryInfo = (hlong*) calloc(NboundaryFaces*(NfaceVertices+1), sizeof(hlong));
  for(hlong n=0;n<Nelements;++n){
    int elementType;
    hlong v1, v2, v3, v4, v5, v6, v7, v8;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);

    if(elementType==3){ // quad boundary face
      sscanf(buf, "%*d%*d %*d" hlongFormat "%*d " hlongFormat hlongFormat hlongFormat hlongFormat,
             boundaryInfo+bcnt*5, &v1, &v2, &v3, &v4);

      boundaryInfo[bcnt*5+1] = v1-1;
      boundaryInfo[bcnt*5+2] = v2-1;
      boundaryInfo[bcnt*5+3] = v3-1;
      boundaryInfo[bcnt*5+4] = v4-1;
      ++bcnt;
    }

    if(elementType==5){  // hex code is 5
      if(start<=Nhexes && Nhexes<=end){
        sscanf(buf,
               "%*d%*d%*d " hlongFormat " %*d"
               hlongFormat hlongFormat hlongFormat hlongFormat hlongFormat hlongFormat hlongFormat hlongFormat,
               elementInfo+cnt,
               &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8);

        EToV[cnt*Nverts+0] = v1-1;
        EToV[cnt*Nverts+1] = v2-1;
        EToV[cnt*Nverts+2] = v3-1;
        EToV[cnt*Nverts+3] = v4-1;
        EToV[cnt*Nverts+4] = v5-1;
        EToV[cnt*Nverts+5] = v6-1;
        EToV[cnt*Nverts+6] = v7-1;
        EToV[cnt*Nverts+7] = v8-1;

        //      printf("%d: %d,%d,%d,%d %d,%d,%d,%d", cnt, v1-1, v2-1,v3-1,v4-1,v5-1,v6-1,v7-1,v8-1);

        ++cnt;
      }
      ++Nhexes;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  NboundaryFaces = bcnt;

  /* record number of found hexes */
  Nelements = (dlong) NhexesLocal;

  /* collect vertices for each element */
  EX = (dfloat*) calloc(Nverts*Nelements, sizeof(dfloat));
  EY = (dfloat*) calloc(Nverts*Nelements, sizeof(dfloat));
  EZ = (dfloat*) calloc(Nverts*Nelements, sizeof(dfloat));
  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<Nverts;++n){
      hlong vid = EToV[e*Nverts+n];
      EX[e*Nverts+n] = VX[vid];
      EY[e*Nverts+n] = VY[vid];
      EZ[e*Nverts+n] = VZ[vid];
    }
  }

  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);
  free(VZ);
}

