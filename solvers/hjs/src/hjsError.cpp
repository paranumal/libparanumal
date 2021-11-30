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
#define CIRCLE_TEST 0
#define SQUARE_TEST 0
#define INTERSECTINGCIRCLE_TEST 1
#define ELLIPSE_TEST 0

#include "hjs.hpp"
void hjs_t::Error(dfloat time, int tstep){

}


// #include "lss.hpp"

// void lss_t::Error(dfloat time, int tstep){
 
// #if 0
//  dlong   Nentries = mesh.Nelements*mesh.Np;
//  const dfloat effective_time = time - 0.1;  //- 100.0*timeStepper->GetTimeStep(); 
 

//  string name;
//  settings.getSetting("OUTPUT FILE NAME", name);
//  char fname[BUFSIZ];
//  sprintf(fname, "%s_%04d.dat", name.c_str(), mesh.N);

//  FILE *fp;

//   fp = fopen(fname, "a");


// dfloat *test_0   = (dfloat *)calloc(Nentries, sizeof(dfloat)); 
// occa::memory o_test = mesh.device.malloc(Nentries*sizeof(dfloat), test_0);

//  o_test.copyFrom(o_phi); 
//  mesh.ogs->GatherScatter(o_test, ogs_dfloat, ogs_add, ogs_sym);
//  linAlg.amx(Nentries, 1.0, o_invDegree, o_test); 
//  o_test.copyTo(phi);

//  // dfloat *linf = (dfloat *)calloc(Nentries, sizeof(dfloat)); 
//  dfloat linf = 0;
//  for(int n=0; n<Nentries; n++ ){
//     const dfloat pn    = phi[n];
//     const dfloat exact = phiex[n];
//     //
//     if( fabs(exact)<effective_time){
//       linf = std::max(linf, fabs(pn-exact)); 
//     }
//   }
  
//   //
//   printf("Writing data\n");
//   fprintf(fp,"%d %.4e %.8e\n", mesh.N, time, linf); 
//   fclose(fp);

// // #endif



//   #else

//    dlong   Nentries = mesh.Nelements*mesh.Np;
//    dfloat *test_0   = (dfloat *)calloc(Nentries, sizeof(dfloat)); 
//    dfloat *test_1   = (dfloat *)calloc(Nentries, sizeof(dfloat)); 
//    dfloat *test_2   = (dfloat *)calloc(Nentries, sizeof(dfloat)); 
//    occa::memory o_test = mesh.device.malloc(Nentries*sizeof(dfloat), test_0);

//    for(int i=0;i<Nentries; i++) test_2[i] = 1.0; 

//   // all ones
//   occa::memory o_test2 = mesh.device.malloc(Nentries*sizeof(dfloat), test_2); 


 

// #if 1
//    o_test.copyFrom(o_phi); 
//    mesh.ogs->GatherScatter(o_test, ogs_dfloat, ogs_add, ogs_sym);
//    linAlg.amx(Nentries, 1.0, o_invDegree, o_test); 
//    o_test.copyTo(phi);
// #endif

//  const dfloat beps0 = 0.10; 
//  const dfloat beps1 = 1.00; 

//   dfloat hmin = mesh.MinCharacteristicLength();

//   //printf("hmin = %.4e\n", hmin);

//  dfloat linfl = 0; dfloat linfg = 0; dfloat Lint = 0; 

//   for(int n=0; n<Nentries; n++){
//     const dfloat xn = mesh.x[n]; 
//     const dfloat yn = mesh.y[n]; 

//     #if CIRCLE_TEST
//     Lint = 2*M_PI*1.0; 
//     const dfloat pn    = phi[n];
//     const dfloat exact = sqrt(xn*xn + yn*yn) - 1.0;
//     // const dfloat init = pow( (xn-1.0)*(xn-1.0) + (yn-1.0)*(yn-1.0) + 0.1, 1.0)*(sqrt(xn*xn + yn*yn) - 1.0) ;
    
//     test_0[n] = 0.0; 
//     test_1[n] = 0.0; 

//     // phi[n] = abs(exact)<= (time+0.05) ? phi[n]:init; 
    
//     if(fabs(exact)<= beps0 ){
//       test_0[n] = fabs(pn-exact);
//       //
//       dfloat tphi = 0.5*(1+tanh(M_PI*pn/hmin)); 
//       dfloat ephi = 0.5*(1+tanh(M_PI*exact/hmin)); 
//       test_2[n]   = tphi-ephi; 


//       linfl     = std::max(linfl, fabs(pn-exact)); 
//     }

//     if(fabs(exact)<= beps1 ){
//       test_1[n] = fabs(pn-exact); 
//       linfg     = std::max(linfg, fabs(pn-exact)); 
//     }
   
//     #elif SQUARE_TEST
//     Lint = 8.0000000; 
//      const dfloat dx = fabs(xn) - 1.0; 
//      const dfloat dy = fabs(yn) - 1.0; 
//      //
//      const dfloat tdx = std::max(dx,0.0); 
//      const dfloat tdy = std::max(dy,0.0); 
//      const dfloat d  = sqrt(tdx*tdx + tdy*tdy); 

//      const dfloat pn    = phi[n];
//      const dfloat exact = d + std::min(std::max(dx,dy),0.0);
    
//     test_0[n] = 0.0; 
//     test_1[n] = 0.0; 
    
   
//     if(fabs(exact)<= beps0 ){
//       test_0[n] = fabs(pn-exact);
//       //
//       dfloat tphi = 0.5*(1+tanh(M_PI*pn/hmin)); 
//       dfloat ephi = 0.5*(1+tanh(M_PI*exact/hmin)); 
//       test_2[n]   = tphi-ephi; 

//       linfl     = std::max(linfl, fabs(pn-exact)); 
//     }

//     if(fabs(exact)<= beps1 ){
//       test_1[n] = fabs(pn-exact); 
//       linfg     = std::max(linfg, fabs(pn-exact)); 
//     }

//     #elif INTERSECTINGCIRCLE_TEST
//  Lint = 9.384775293622596; 
//     const dfloat pn    = phi[n];
//     const dfloat exact = phiex[n];
    
//     test_0[n] = 0.0; 
//     test_1[n] = 0.0; 
    
//     if(fabs(exact)<= beps0 ){
//       test_0[n] = fabs(pn-exact);
//       //
//       dfloat tphi = 0.5*(1+tanh(M_PI*pn/hmin)); 
//       dfloat ephi = 0.5*(1+tanh(M_PI*exact/hmin)); 
//       test_2[n]   = tphi-ephi; 

//       linfl     = std::max(linfl, fabs(pn-exact)); 
//     }

//     if(fabs(exact)<= beps1 ){
//       test_1[n] = fabs(pn-exact); 
//       linfg     = std::max(linfg, fabs(pn-exact)); 
//     }




//     #elif ELLIPSE_TEST
//     Lint = 4.844224110273838099; 
//     const dfloat pn    = phi[n];
//     const dfloat exact = phiex[n];
    
//     test_0[n] = 0.0; 
//     test_1[n] = 0.0; 
    
//     if(fabs(exact)<= beps0 ){
//       test_0[n] = fabs(pn-exact);
//       //
//       dfloat tphi = 0.5*(1+tanh(M_PI*pn/hmin)); 
//       dfloat ephi = 0.5*(1+tanh(M_PI*exact/hmin)); 
//       test_2[n]   = tphi-ephi; 

//       linfl     = std::max(linfl, fabs(pn-exact)); 
//     }

//     if(fabs(exact)<= beps1 ){
//       test_1[n] = fabs(pn-exact); 
//       linfg     = std::max(linfg, fabs(pn-exact)); 
//     }
//     #endif 
//   }



//   o_test.copyFrom(test_1); 
//   MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_test, o_Mq);
//   dfloat normg2 = sqrt(linAlg.innerProd(Nentries, o_test, o_Mq, comm));

//   o_test.copyFrom(test_0); 
//   MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_test, o_Mq);
//   dfloat norml2 = sqrt(linAlg.innerProd(Nentries, o_test, o_Mq, comm));

//   o_test.copyFrom(test_0); 
//   MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_test, o_Mq);
//   dfloat norml1 = linAlg.innerProd(Nentries, o_test2, o_Mq, comm)/Lint;



// if(mesh.rank==0){
//     printf("%5.2f (%d), %.4e %.4e %.4e %.4e %.4e(time, timestep, L2_g, L2_l,Linf_g, Linf_l, L1)\n", 
//       time, tstep, normg2, norml2, linfg,linfl, norml1);  

// #if 0
//     string name;
//     settings.getSetting("OUTPUT FILE NAME", name);
//     char fname[BUFSIZ];
//     sprintf(fname, "ellipse_h2_N_%d_eps_%s.dat", mesh.N, name.c_str());

//    FILE *fp;
//    fp = fopen(fname, "a");
  
//    fprintf(fp, "%.4e %d %.4e %.4e %.4e %.4e %.4e\n", time, tstep, normg2, norml2, linfg,linfl, norml1);  
//     fclose(fp);
// #endif


//    }



//    // for(int i=0; i<Nentries; i++){
//    //  phi[i] = abs(phi[i])<=time ? phi[i]:100; 
//    // }


//    // if(mesh.rank==0){
//    //  printf("%5.2f (%d), %.8e %.8e (time, timestep, norm_g, norm_l)\n", time, tstep, normg2, norml2);    
//    // }

//    free(test_0); free(test_1);


//   // if(mesh.rank==0)
//   //   printf("%5.2f (%d), %.8e (time, timestep, norm)\n", time, tstep, norm2);

//   // if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

//   //   o_q.copyTo(q);
//   //   // o_phi.copyTo(phi); 
//   //   if(redistance){
//   //     // o_sgnq.copyTo(sgnq);
//   //     subcell->o_ElementList.copyTo(subcell->ElementList);
//   //   }   
//   //   // output field files
//   //   string name;
//   //   settings.getSetting("OUTPUT FILE NAME", name);
//   //   char fname[BUFSIZ];
//   //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

//   //   PlotFields(q, fname);
//   // }
// #endif
// }
