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

#include "elliptic.h"

void ellipticScaledAdd(elliptic_t *elliptic, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  mesh_t *mesh = elliptic->mesh;

  dlong Ntotal = mesh->Nelements*mesh->Np;

  // b[n] = alpha*a[n] + beta*b[n] n\in [0,Ntotal)
  occaTimerTic(mesh->device,"scaledAddKernel");
  elliptic->scaledAddKernel(Ntotal, alpha, o_a, beta, o_b);
  occaTimerToc(mesh->device,"scaledAddKernel");
}

dfloat ellipticWeightedInnerProduct(elliptic_t *elliptic, occa::memory &o_w, occa::memory &o_a, occa::memory &o_b){

  mesh_t *mesh = elliptic->mesh;
  dfloat *tmp = elliptic->tmp;
  dlong Nblock = elliptic->Nblock;
  dlong Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = elliptic->o_tmp;

  occaTimerTic(mesh->device,"weighted inner product2");
  if(elliptic->options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    elliptic->weightedInnerProduct2Kernel(Ntotal, o_w, o_a, o_b, o_tmp);
  else
    elliptic->innerProductKernel(Ntotal, o_a, o_b, o_tmp);

  occaTimerToc(mesh->device,"weighted inner product2");

  o_tmp.copyTo(tmp);

  dfloat wab = 0;
  for(dlong n=0;n<Nblock;++n){
    wab += tmp[n];
  }

  dfloat globalwab = 0;
  MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

  return globalwab;
}

dfloat ellipticWeightedNorm2(elliptic_t *elliptic, occa::memory &o_w, occa::memory &o_a){

  mesh_t *mesh = elliptic->mesh;
  dfloat *tmp = elliptic->tmp;
  dlong Nblock = elliptic->Nblock;
  dlong Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = elliptic->o_tmp;

  occaTimerTic(mesh->device,"weighted inner product2");
  if(elliptic->options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    elliptic->weightedNorm2Kernel(Ntotal, o_w, o_a, o_tmp);
  else
    elliptic->norm2Kernel(Ntotal, o_a, o_tmp);

  occaTimerToc(mesh->device,"weighted norm2");

  o_tmp.copyTo(tmp);

  dfloat wab = 0;
  for(dlong n=0;n<Nblock;++n){
    wab += tmp[n];
  }

  dfloat globalwab = 0;
  MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

  return globalwab;
}


dfloat ellipticInnerProduct(elliptic_t *elliptic, occa::memory &o_a, occa::memory &o_b){

  mesh_t *mesh = elliptic->mesh;
  dfloat *tmp = elliptic->tmp;
  dlong Nblock = elliptic->Nblock;
  dlong Ntotal = mesh->Nelements*mesh->Np;

  occa::memory &o_tmp = elliptic->o_tmp;

  occaTimerTic(mesh->device,"inner product");
  elliptic->innerProductKernel(Ntotal, o_a, o_b, o_tmp);
  occaTimerToc(mesh->device,"inner product");

  o_tmp.copyTo(tmp);

  dfloat ab = 0;
  for(dlong n=0;n<Nblock;++n){
    ab += tmp[n];
  }

  dfloat globalab = 0;
  MPI_Allreduce(&ab, &globalab, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

  return globalab;
}
