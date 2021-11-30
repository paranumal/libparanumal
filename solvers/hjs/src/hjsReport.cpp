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
#define CIRCLE_TEST 1
#define SQUARE_TEST 0
#define INTERSECTINGCIRCLE_TEST 0

#include "hjs.hpp"

void hjs_t::Report(dfloat time, int tstep){

static int frame=0;

  //compute q.M*q
  mesh.MassMatrixApply(o_q, o_Mq);

  dlong Nentries = mesh.Nelements*mesh.Np;
  dfloat norm2 = sqrt(platform.linAlg.innerProd(Nentries, o_q, o_Mq, mesh.comm));

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", time, tstep, norm2);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // copy data back to host
    o_q.copyTo(q);

    // output field files
    string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

    PlotFields(q, fname);
  }

	
}



// #include "lss.hpp"

// void lss_t::Report(dfloat time, int tstep){

//   Error(time, tstep);

//   static int frame=0;
// // #if 0
// //   // //compute q.M*q
// //   // lss version
// //   MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_phi, o_Mq);

// //   dlong Nentries = mesh.Nelements*mesh.Np;
// //   dfloat norm2 = sqrt(linAlg.innerProd(Nentries, o_phi, o_Mq, comm));
// // #else
// //    dlong Nentries = mesh.Nelements*mesh.Np;
// //    dfloat *test = (dfloat *)calloc(Nentries, sizeof(dfloat)); 

// //    occa::memory o_test = mesh.device.malloc(Nentries*sizeof(dfloat), test);
   
// //    o_test.copyFrom(o_phi); 
// //    mesh.ogs->GatherScatter(o_test, ogs_dfloat, ogs_add, ogs_sym);
// //    linAlg.amx(Nentries, 1.0, o_invDegree, o_test); 

// //    o_test.copyTo(phi);

// //   for(int n=0; n<Nentries; n++){
// //     const dfloat xn = mesh.x[n]; 
// //     const dfloat yn = mesh.y[n]; 

// //     #if CIRCLE_TEST
// //         dfloat exact = sqrt(xn*xn + yn*yn) - 1.0;
// //         test[n] = fabs(phi[n] - exact); 
// //     #elif SQUARE_TEST


// //     #elif INTERSECTINGCIRCLE_TEST



// //     #endif 
// //   }

// //   o_test.copyFrom(test); 

// //   MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_test, o_Mq);
// //   dfloat norm2 = sqrt(linAlg.innerProd(Nentries, o_test, o_Mq, comm));

// //   // norm2 *= norm2; // area 
  
// //   // error for r=1 circle
// //   // norm2 = M_PI - norm2;







// //    free(test); 

// // #endif

//   // if(mesh.rank==0)
//   //   printf("%5.2f (%d), %.8e (time, timestep, norm)\n", time, tstep, norm2);

//   if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

//     o_q.copyTo(q);
//     // o_phi.copyTo(phi); 
//     if(redistance){
//       // o_sgnq.copyTo(sgnq);
//       subcell->o_ElementList.copyTo(subcell->ElementList);
//     }   
//     // output field files
//     string name;
//     settings.getSetting("OUTPUT FILE NAME", name);
//     char fname[BUFSIZ];
//     sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

//     PlotFields(q, fname);
//   }
// }
