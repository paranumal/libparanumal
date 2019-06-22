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

#include "ins.h"
// #include "ogsInterface.h"

dfloat insMean(ins_t *ins, occa::memory o_q){

  mesh_t *mesh = ins->mesh;
  dfloat qmeanLocal;
  dfloat qmeanGlobal;
  
 // dlong Nblock = (mesh->cubNp*mesh->Nelements+blockSize-1)/blockSize;
  dlong Nblock = ins->Nblock; 
  dfloat *tmp = (dfloat *)calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  occa::memory o_tmp = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), tmp);

  mesh->sumKernel(mesh->Nelements*mesh->Np, o_q, o_tmp);
  o_tmp.copyTo(tmp);
  // finish reduction
  qmeanLocal = 0;
  for(dlong n=0;n<Nblock;++n)
    qmeanLocal += tmp[n];

  // globalize reduction
  MPI_Allreduce(&qmeanLocal, &qmeanGlobal, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
  qmeanGlobal /= ((dfloat) ins->pSolver->NelementsGlobal*(dfloat)mesh->Np);

  printf("Mean of the field : %.8e\n", qmeanGlobal); 

  return qmeanGlobal; 
}
