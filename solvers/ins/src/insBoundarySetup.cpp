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

#include "ins.hpp"

void ins_t::BoundarySetup(){

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  mapB = (int *) calloc(mesh.Nelements*mesh.Np,sizeof(int));
  const int largeNumber = 1<<20;
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int n=0;n<mesh.Np;n++) mapB[n+e*mesh.Np] = largeNumber;
    for (int f=0;f<mesh.Nfaces;f++) {
      int bc = mesh.EToB[f+e*mesh.Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh.Nfp;n++) {
          int fid = mesh.faceNodes[n+f*mesh.Nfp];
          mapB[fid+e*mesh.Np] = mymin(bc,mapB[fid+e*mesh.Np]);
        }
      }
    }
  }
  mesh.ogs->GatherScatter(mapB, ogs_int, ogs_min, ogs_sym);

  for (dlong n=0;n<mesh.Nelements*mesh.Np;n++) {
    if (mapB[n] == largeNumber) {//no boundary
      mapB[n] = 0.;
    }
  }

  o_mapB = mesh.device.malloc(mesh.Nelements*mesh.Np*sizeof(int), mapB);
}
