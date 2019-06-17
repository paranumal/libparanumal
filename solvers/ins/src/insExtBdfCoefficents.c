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

void insExtBdfCoefficents(ins_t *ins, int order) {

  if(order==1) {
     //advection, first order in time, increment
    ins->g0 =  1.0f; 
    dfloat extbdfB[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfA[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};
    
    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder = 1;    
    
    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0; 
  } else if(order==2) {
    //advection, second order in time, increment
    ins->g0 =  1.5f;
    dfloat extbdfB[3] = {2.0f,-0.5f, 0.0f};
    dfloat extbdfA[3] = {2.0f,-1.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};

    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder=2;

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0; 
  } else if(order==3) {
    //advection, third order in time, increment
    ins->g0 =  11.f/6.f;
    dfloat extbdfB[3] = {3.0f,-1.5f, 1.0f/3.0f};
    dfloat extbdfA[3] = {3.0f,-3.0f, 1.0f};
    dfloat extbdfC[3] = {2.0f,-1.0f, 0.0f};
    
    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder=3;

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0; 
  }
}
