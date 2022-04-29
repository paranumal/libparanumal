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

  nonPmlElements.malloc(NnonPmlElements);
  pmlElements.malloc(NpmlElements);
  pmlIds.malloc(NpmlElements);

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

  o_pmlElements = platform.malloc<dlong>(pmlElements);
  o_pmlIds = platform.malloc<dlong>(pmlIds);
  o_nonPmlElements = platform.malloc<dlong>(nonPmlElements);
}


void mesh_t::MultiRatePmlSetup(){

  mrNnonPmlElements.malloc(mrNlevels, 0);
  mrNpmlElements.malloc(mrNlevels, 0);

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

  mrNonPmlElements.malloc(mrNlevels);
  mrPmlElements.malloc(mrNlevels);
  mrPmlIds.malloc(mrNlevels);
  for (int lev=0;lev<mrNlevels;lev++) {
    mrNonPmlElements[lev].malloc(mrNnonPmlElements[lev]);
    mrPmlElements[lev].malloc(mrNpmlElements[lev]);
    mrPmlIds[lev].malloc(mrNpmlElements[lev]);

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

  o_mrNonPmlElements.malloc(mrNlevels);
  o_mrPmlElements.malloc(mrNlevels);
  o_mrPmlIds.malloc(mrNlevels);

  for (int lev=0;lev<mrNlevels;lev++){
    o_mrPmlElements[lev]   = platform.malloc<dlong>(mrPmlElements[lev]);
    o_mrPmlIds[lev] = platform.malloc<dlong>(mrPmlIds[lev]);
    o_mrNonPmlElements[lev] = platform.malloc<dlong>(mrNonPmlElements[lev]);
  }
}

} //namespace libp
