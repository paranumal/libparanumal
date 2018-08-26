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

void insVelocityUpdate(ins_t *ins, dfloat time, int stage, 
                        occa::memory o_rkGP,
                        occa::memory o_rkU){

  mesh_t *mesh = ins->mesh;

  // U^s = Uhat - dt * GP^s + dt*\sum^s-1 pa_si GP^i
  occaTimerTic(mesh->device,"VelocityUpdate");
  ins->velocityUpdateKernel(mesh->Nelements,
                              stage,
                              ins->ARKswitch,
                              ins->dt,
                              ins->fieldOffset,
                              ins->o_prkA,
                              ins->o_prkB,
                              o_rkGP,
                              ins->o_GP,
                              o_rkU);
  occaTimerToc(mesh->device,"VelocityUpdate");
}
