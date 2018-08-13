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

#include "mppf.h"

void mppfPressureRhs(mppf_t *mppf, dfloat time, occa::memory o_rkU){

  mesh_t *mesh = mppf->mesh;


  // Compute Explicit Diffusive Terms; i.e. DU = curl*curl*u^*  and  SU = grad(mu)\cdot D(u^*)
  mppfExplicitDiffusive(mppf, time, mppf->o_U, mppf->o_DU, mppf->o_SU); 


// Give exact DU and SU are fine
#if 0
  for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      dfloat dux = 2*M_PI*M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(time);
      dfloat duy =-2*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(time);

      mppf->rkU[id + 0*mppf->fieldOffset] = dux;
      mppf->rkU[id + 1*mppf->fieldOffset] = duy;
    }
  }
  
  mppf->o_DU.copyFrom(mppf->rkU); 


  for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      dfloat phix   = -M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(time);
      dfloat phiy   = -M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(time);

      dfloat ux =  M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(time);
      dfloat uy = -M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(time);
      dfloat vx =  M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(time);
      dfloat vy = -M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(time);

      dfloat mux = 0.5*(mppf->mu1-mppf->mu2)*phix;
      dfloat muy = 0.5*(mppf->mu1-mppf->mu2)*phiy;


      dfloat sux = mux*(2.0*ux) + muy*(uy+ vx);
      dfloat suy = mux*(vx + uy) + muy*(2.0*vy);

      mppf->rkU[id + 0*mppf->fieldOffset] = sux;
      mppf->rkU[id + 1*mppf->fieldOffset] = suy;
    }
  }

  mppf->o_SU.copyFrom(mppf->rkU); 

#endif

  // Upadte i.e. o_rkU = -NU + U^/dt + (1/rh0 - 1/rho^n+1)*GP^(*,n+1) -mu^(n+1)/rho^(n+1)*DU 
  //                     + 1/rho^(n+1)*SU - lamdba*Psi*Ghi + 1/(rho^n+1)*f^(n+1) 
  mppf->explicitUpdateKernel(mesh->Nelements,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mppf->dt,
                           mppf->time,
                           mppf->chA,
                           mppf->o_extbdfA,
                           mppf->o_extbdfB,
                           mppf->fieldOffset,
                           mppf->o_U,
                           mppf->o_Rho,
                           mppf->o_Mu,
                           mppf->o_NU,
                           mppf->o_DU,
                           mppf->o_SU,
                           mppf->o_GP,
                           mppf->o_Psi,
                           mppf->o_Phi,
                           mppf->o_GPhi,
                           o_rkU);


// Give exact NU
#if 0 // Divergence is fine
  for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      dfloat rhsX = cos(time)*(cos(M_PI*y)*sin(M_PI*x) + M_PI*cos(M_PI*x)*sin(M_PI*y));
      dfloat rhsY = -cos(time)*(cos(M_PI*x)*sin(M_PI*y) - M_PI*cos(M_PI*y)*sin(M_PI*x));
 

     
      mppf->rkU[id + 0*mppf->fieldOffset] = mppf->dt*rhsX;
      mppf->rkU[id + 1*mppf->fieldOffset] = mppf->dt*rhsY;
      // mppf->rhsP[id] = rhsp;
    }
  }

  o_rkU.copyFrom(mppf->rkU); 

#endif

  
  // rhsP = Div Uhat
  mppfDivergence(mppf, time, o_rkU, mppf->o_rhsP);



// Give exact NU
#if 0
  for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      dfloat rhsp = -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(time);
 

     
      mppf->rhsP[id] = mppf->dt*rhsp;
      // mppf->rhsP[id] = rhsp;
    }
  }

  mppf->o_rhsP.copyFrom(mppf->rhsP); 

#endif




  
  // // rhsP = -MM*rho0* Div Uhat
  occaTimerTic(mesh->device,"PoissonRhsForcing");
  mppf->pressureRhsKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mppf->idt,
                              mesh->o_MM, 
                              mppf->o_rhsP);
  occaTimerToc(mesh->device,"PoissonRhsForcing");

  
  

#if 0
  //add penalty from jumps in previous pressure
  mppf->poissonPenaltyKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_vgeo,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                mppf->tau,
                                mesh->o_x,
                                mesh->o_y,
                                t,
                                mppf->dt,
                                mppf->c0,
                                mppf->c1,
                                mppf->c2,
                                mppf->index,
                                (mesh->Nelements+mesh->totalHaloPairs),
                                mppf->o_P,
                                mppf->o_rhsP);
  #endif







}
