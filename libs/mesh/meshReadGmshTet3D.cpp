/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

namespace libp {

/*
   purpose: read gmsh tetrahedra mesh
*/
void mesh_t::ReadGmshTet3D(const std::string fileName){

  FILE *fp = fopen(fileName.c_str(), "r");
  LIBP_ABORT("Cannot open file: " << fileName,
             fp==NULL);

  char buf[BUFSIZ];
  do{
    //read to end of line
    LIBP_ABORT("Error reading mesh file: " << fileName,
               !fgets(buf, BUFSIZ, fp));
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  //read to end of line
  LIBP_ABORT("Error reading mesh file: " << fileName,
             !fgets(buf, BUFSIZ, fp));
  sscanf(buf, hlongFormat, &(Nnodes));

  /* allocate space for node coordinates */
  memory<dfloat> VX(Nnodes);
  memory<dfloat> VY(Nnodes);
  memory<dfloat> VZ(Nnodes);

  /* load nodes */
  for(hlong n=0;n<Nnodes;++n){
    //read to end of line
    LIBP_ABORT("Error reading mesh file: " << fileName,
               !fgets(buf, BUFSIZ, fp));
    sscanf(buf, "%*d" dfloatFormat dfloatFormat dfloatFormat,
           VX.ptr()+n, VY.ptr()+n, VZ.ptr()+n);
  }

  /* look for section with Element node data */
  do{
    //read to end of line
    LIBP_ABORT("Error reading mesh file: " << fileName,
               !fgets(buf, BUFSIZ, fp));
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  hlong gNelements;
  //read to end of line
  LIBP_ABORT("Error reading mesh file: " << fileName,
             !fgets(buf, BUFSIZ, fp));
  sscanf(buf, hlongFormat, &gNelements);

  /* find # of tets */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  hlong Ntets = 0, gNboundaryFaces = 0;
  for(hlong n=0;n<gNelements;++n){
    int ElementType;
    //read to end of line
    LIBP_ABORT("Error reading mesh file: " << fileName,
               !fgets(buf, BUFSIZ, fp));
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
  hlong start = rank*chunk + std::min(rank, remainder);
  hlong end = start + NtetsLocal-1;

  /* allocate space for Element node index data */
  EToV.malloc(NtetsLocal*Nverts);
  elementInfo.malloc(NtetsLocal);

  /* scan through file looking for tetrahedra elements */
  hlong cnt=0, bcnt = 0;
  Ntets = 0;

  boundaryInfo.malloc(gNboundaryFaces*4);
  for(hlong n=0;n<gNelements;++n){
    int ElementType;
    hlong v1, v2, v3, v4;
    //read to end of line
    LIBP_ABORT("Error reading mesh file: " << fileName,
               !fgets(buf, BUFSIZ, fp));
    sscanf(buf, "%*d%d", &ElementType);
    if(ElementType==2){ // boundary face
      sscanf(buf, "%*d%*d %*d" hlongFormat "%*d" hlongFormat hlongFormat hlongFormat,
             boundaryInfo.ptr()+bcnt*4, &v1, &v2, &v3);
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
               elementInfo.ptr()+cnt,&v1, &v2, &v3, &v4);
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
  EX.malloc(Nverts*Nelements);
  EY.malloc(Nverts*Nelements);
  EZ.malloc(Nverts*Nelements);
  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<Nverts;++n){
      hlong vid = EToV[e*Nverts+n];
      EX[e*Nverts+n] = VX[vid];
      EY[e*Nverts+n] = VY[vid];
      EZ[e*Nverts+n] = VZ[vid];
    }
  }
}

} //namespace libp
