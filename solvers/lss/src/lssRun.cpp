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

#include "lss.hpp"

void lss_t::Run(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  initialConditionKernel(mesh.Nelements,
                         startTime,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_q);
  if(advection){
  setFlowFieldKernel(mesh.Nelements,
                    offset, 
                    startTime,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_U, 
                    o_q);
}else{
  regularizedSignKernel(mesh.Nelements, 
                       startTime,
                       eps, 
                       mesh.o_x,
                       mesh.o_y,
                       mesh.o_z,
                       o_q,
                       o_sgnq); 

 if(subcellStabilization){


   const int all = 1; // all elements, do not use element list
    projectKernel(mesh.Nelements,
                  all,  
                  subcell->o_ElementList, 
                  subcell->o_PMT,
                  o_sgnq, 
                  o_ssgnq);


    // skylineKernel(mesh.Nelements,
    //             mesh.o_ggeo,
    //             subcell->o_ModMap, 
    //             mesh.o_MM,
    //             subcell->o_invVT,
    //             subcell->o_LSF,
    //             subcell->o_BLD, 
    //             o_q,
    //             subcell->o_ElementList); 

    skyline1DKernel(mesh.Nelements,
                mesh.o_sgeo,
                subcell->o_edgeNodes, 
                subcell->o_MM1D,
                subcell->o_invVT1D,
                subcell->o_LSF,
                subcell->o_BLD, 
                o_q,
                subcell->o_ElementList); 
  
  }

}

 // dfloat outputInterval; 
 // settings.getSetting("OUTPUT INTERVAL", outputInterval);
 finalTime = startTime + 1*timeStepper->GetTimeStep();
 // finalTime = startTime + 53*outputInterval;
  timeStepper->Run(o_q, startTime, finalTime);
}
