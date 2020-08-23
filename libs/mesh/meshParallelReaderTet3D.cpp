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
#include "mesh3D.hpp"

/*
   purpose: read gmsh tetrahedra mesh
*/
void meshTet3D::ParallelReader(const char *fileName){

  FILE *fp = fopen(fileName, "r");

  dim = 3;
  Nverts = 4; // number of vertices per element
  Nfaces = 4;

  // vertices on each face
  int faceVertices_[4][3] = {{0,1,2},{0,1,3},{1,2,3},{2,0,3}};

  NfaceVertices = 3;
  faceVertices =
    (int*) calloc(NfaceVertices*Nfaces, sizeof(int));
  memcpy(faceVertices, faceVertices_[0], 12*sizeof(int));

  if(fp==NULL){
    stringstream ss;
    ss << "Cannot open file: " << fileName;
    LIBP_ABORT(ss.str())
  }

  char buf[BUFSIZ];
  do{
    if (!fgets(buf, BUFSIZ, fp)) { //read to end of line
      stringstream ss;
      ss << "Error reading mesh file: " << fileName;
      LIBP_ABORT(ss.str())
    }
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  if (!fgets(buf, BUFSIZ, fp)) { //read to end of line
    stringstream ss;
    ss << "Error reading mesh file: " << fileName;
    LIBP_ABORT(ss.str())
  }
  sscanf(buf, hlongFormat, &(Nnodes));

  /* allocate space for node coordinates */
  dfloat *VX = (dfloat*) calloc(Nnodes, sizeof(dfloat));
  dfloat *VY = (dfloat*) calloc(Nnodes, sizeof(dfloat));
  dfloat *VZ = (dfloat*) calloc(Nnodes, sizeof(dfloat));

  /* load nodes */
  for(hlong n=0;n<Nnodes;++n){
    if (!fgets(buf, BUFSIZ, fp)) { //read to end of line
      stringstream ss;
      ss << "Error reading mesh file: " << fileName;
      LIBP_ABORT(ss.str())
    }
    sscanf(buf, "%*d" dfloatFormat dfloatFormat dfloatFormat,
           VX+n, VY+n, VZ+n);
  }

  /* look for section with Element node data */
  do{
    if (!fgets(buf, BUFSIZ, fp)) { //read to end of line
      stringstream ss;
      ss << "Error reading mesh file: " << fileName;
      LIBP_ABORT(ss.str())
    }
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  hlong gNelements;
  if (!fgets(buf, BUFSIZ, fp)) { //read to end of line
    stringstream ss;
    ss << "Error reading mesh file: " << fileName;
    LIBP_ABORT(ss.str())
  }
  sscanf(buf, hlongFormat, &gNelements);

  /* find # of tets */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  hlong Ntets = 0, gNboundaryFaces = 0;
  for(hlong n=0;n<gNelements;++n){
    int ElementType;
    if (!fgets(buf, BUFSIZ, fp)) { //read to end of line
      stringstream ss;
      ss << "Error reading mesh file: " << fileName;
      LIBP_ABORT(ss.str())
    }
    sscanf(buf, "%*d%d", &ElementType);
    if(ElementType==4) ++Ntets; // tet code is 4
    if(ElementType==2) ++gNboundaryFaces;
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  hlong chunk = (hlong) Ntets/size;
  int remainder = (int) (Ntets - chunk*size);

  hlong NtetsLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  hlong start = rank*chunk + mymin(rank, remainder);
  hlong end = start + NtetsLocal-1;

  /* allocate space for Element node index data */

  EToV
    = (hlong*) calloc(NtetsLocal*Nverts, sizeof(hlong));
  elementInfo
    = (hlong*) calloc(NtetsLocal,sizeof(hlong));

  /* scan through file looking for tetrahedra elements */
  hlong cnt=0, bcnt = 0;
  Ntets = 0;

  boundaryInfo = (hlong*) calloc(gNboundaryFaces*4, sizeof(hlong));
  for(hlong n=0;n<gNelements;++n){
    int ElementType;
    hlong v1, v2, v3, v4;
    if (!fgets(buf, BUFSIZ, fp)) { //read to end of line
      stringstream ss;
      ss << "Error reading mesh file: " << fileName;
      LIBP_ABORT(ss.str())
    }
    sscanf(buf, "%*d%d", &ElementType);
    if(ElementType==2){ // boundary face
      sscanf(buf, "%*d%*d %*d" hlongFormat "%*d" hlongFormat hlongFormat hlongFormat,
             boundaryInfo+bcnt*4, &v1, &v2, &v3);
      boundaryInfo[bcnt*4+1] = v1-1;
      boundaryInfo[bcnt*4+2] = v2-1;
      boundaryInfo[bcnt*4+3] = v3-1;
      ++bcnt;
    }

    if(ElementType==4){  // tet code is 4
      if(start<=Ntets && Ntets<=end){
        sscanf(buf,
               "%*d%*d%*d " hlongFormat " %*d"
               hlongFormat hlongFormat hlongFormat hlongFormat,
               elementInfo+cnt,&v1, &v2, &v3, &v4);
        /* read vertex triplet for trianngle */
        EToV[cnt*Nverts+0] = v1-1;
        EToV[cnt*Nverts+1] = v2-1;
        EToV[cnt*Nverts+2] = v3-1;
        EToV[cnt*Nverts+3] = v4-1;
        ++cnt;
      }
      ++Ntets;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  NboundaryFaces = bcnt;

  /* record number of found tets */
  Nelements = (dlong) NtetsLocal;

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

