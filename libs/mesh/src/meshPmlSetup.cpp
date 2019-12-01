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

void mesh_t::PmlSetup(){

  NnonPmlElements=0;
  NpmlElements=0;

  //count PML elements
  for (dlong e=0;e<Nelements;e++) {
    hlong type = elementInfo[e];

    //these element info flags are reserved for PML elements of certain alignment
    // 100 - x   PML
    // 200 - y   PML
    // 300 - xy  PML
    // 400 - z   PML
    // 500 - xz  PML
    // 600 - yz  PML
    // 700 - xyz PML
    if ((type==100)||(type==200)||(type==300)||
        (type==400)||(type==500)||(type==600)||
        (type==700) )
      NpmlElements++;
    else
      NnonPmlElements++;
  }

  nonPmlElements = (dlong *) malloc(NnonPmlElements*sizeof(dlong));
  pmlElements    = (dlong *) malloc(NpmlElements*sizeof(dlong));
  pmlIds         = (dlong *) malloc(NpmlElements*sizeof(dlong*));

  NnonPmlElements=0;
  NpmlElements=0;
  dlong pmlCnt=0;
  for (dlong e=0;e<Nelements;e++) {
    hlong type = elementInfo[e];

    if ((type==100)||(type==200)||(type==300)||
        (type==400)||(type==500)||(type==600)||
        (type==700) ) {
      pmlElements[NpmlElements] = e;
      pmlIds[NpmlElements++] = pmlCnt++;
    } else
      nonPmlElements[NnonPmlElements++] = e;
  }

  if (NpmlElements) {
    o_pmlElements = device.malloc(NpmlElements*sizeof(dlong), pmlElements);
    o_pmlIds = device.malloc(NpmlElements*sizeof(dlong), pmlIds);
  }

  if (NnonPmlElements)
    o_nonPmlElements = device.malloc(NnonPmlElements*sizeof(dlong), nonPmlElements);
}


void mesh_t::MultiRatePmlSetup(){

  mrNnonPmlElements = (dlong *) calloc(mrNlevels,sizeof(dlong));
  mrNpmlElements    = (dlong *) calloc(mrNlevels,sizeof(dlong));

  //count PML elements
  for (dlong e=0;e<Nelements;e++) {
    hlong type = elementInfo[e];
    int lev = mrLevel[e];

    //these element info flags are reserved for PML elements of certain alignment
    // 100 - x   PML
    // 200 - y   PML
    // 300 - xy  PML
    // 400 - z   PML
    // 500 - xz  PML
    // 600 - yz  PML
    // 700 - xyz PML
    if ((type==100)||(type==200)||(type==300)||
        (type==400)||(type==500)||(type==600)||
        (type==700) )
      for (int l=lev;l<mrNlevels;l++) mrNpmlElements[l]++;
    else
      for (int l=lev;l<mrNlevels;l++) mrNnonPmlElements[l]++;
  }

  mrNonPmlElements = (dlong **) malloc(mrNlevels*sizeof(dlong*));
  mrPmlElements    = (dlong **) malloc(mrNlevels*sizeof(dlong*));
  mrPmlIds         = (dlong **) malloc(mrNlevels*sizeof(dlong*));
  for (int lev=0;lev<mrNlevels;lev++) {
    mrNonPmlElements[lev] = (dlong *) malloc(mrNnonPmlElements[lev]*sizeof(dlong));
    mrPmlElements[lev]    = (dlong *) malloc(mrNpmlElements[lev]*sizeof(dlong));
    mrPmlIds[lev]         = (dlong *) malloc(mrNpmlElements[lev]*sizeof(dlong));

    //reset
    mrNpmlElements[lev] = 0;
    mrNnonPmlElements[lev] = 0;
  }

  dlong pmlCnt=0;
  for (dlong e=0;e<Nelements;e++) {
    hlong type = elementInfo[e];
    int lev = mrLevel[e];

    if ((type==100)||(type==200)||(type==300)||
        (type==400)||(type==500)||(type==600)||
        (type==700) ) {
      for (int l=lev;l<mrNlevels;l++) {
        mrPmlElements[l][mrNpmlElements[l]] = e;
        mrPmlIds[l][mrNpmlElements[l]++] = pmlCnt;
      }
      pmlCnt++;
    } else
      for (int l=lev;l<mrNlevels;l++)
        mrNonPmlElements[l][mrNnonPmlElements[l]++] = e;
  }

  o_mrNonPmlElements = new occa::memory[mrNlevels];
  o_mrPmlElements    = new occa::memory[mrNlevels];
  o_mrPmlIds         = new occa::memory[mrNlevels];

  for (int lev=0;lev<mrNlevels;lev++){
    if (mrNpmlElements[lev]) {
      o_mrPmlElements[lev]   = device.malloc(mrNpmlElements[lev]*sizeof(dlong), mrPmlElements[lev]);
      o_mrPmlIds[lev] = device.malloc(mrNpmlElements[lev]*sizeof(dlong), mrPmlIds[lev]);
    }
    if (mrNnonPmlElements[lev])
      o_mrNonPmlElements[lev] = device.malloc(mrNnonPmlElements[lev]*sizeof(dlong), mrNonPmlElements[lev]);
  }
}