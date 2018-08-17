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

      dfloat dux =  2.0*M_PI*M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(time);
      dfloat duy = -2.0*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(time);

      dfloat rho = 2.0 - cos(M_PI*x)*cos(M_PI*y)*sin(time);
      dfloat mu  = 3.0/200.0 - (cos(M_PI*x)*cos(M_PI*y)*sin(time))/200.0;

      dfloat sux =  (M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*x)*sin(time)*sin(time))/100.0;
      dfloat suy = -(M_PI*M_PI*cos(M_PI*x)*cos(M_PI*x)*cos(M_PI*y)*sin(M_PI*y)*sin(time)*sin(time))/100.0;

      dfloat px = M_PI*cos(M_PI*x)*sin(M_PI*y)*cos(time);
      dfloat py = M_PI*cos(M_PI*y)*sin(M_PI*x)*cos(time);
      
      dfloat phi     = cos(M_PI*x)*cos(M_PI*y)*sin(time);
      dfloat lapPhi = -2.0*M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(time) + mppf->chA*phi;
      dfloat phix = -M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(time);
      dfloat phiy = -M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(time);

      dfloat nux = -(M_PI*sin(2.0*M_PI*x)*(cos(2.0*time)/2.0 - 0.5))/2.0;
      dfloat nuy = -(M_PI*sin(2.0*M_PI*y)*(cos(2.0*time)/2.0 - 0.5))/2.0;

      mppf->DU[id + 0*mppf->fieldOffset]   = dux;
      mppf->DU[id + 1*mppf->fieldOffset]   = duy;

      mppf->SU[id + 0*mppf->fieldOffset]   = sux;
      mppf->SU[id + 1*mppf->fieldOffset]   = suy;

      mppf->Rho[id                     ]   = rho;
      mppf->Mu[id                      ]   = mu;

      mppf->GP[id + 0*mppf->fieldOffset]   = px;
      mppf->GP[id + 1*mppf->fieldOffset]   = py;

      mppf->Psi[id                     ]   = lapPhi;

      mppf->GPhi[id + 0*mppf->fieldOffset] = phix;
      mppf->GPhi[id + 1*mppf->fieldOffset] = phiy;

      mppf->NU[id + 0*mppf->fieldOffset] = nux;
      mppf->NU[id + 1*mppf->fieldOffset] = nuy;



      // mppf->rkU[id + 0*mppf->fieldOffset] = -nux + (1.0/mppf->rho0 - 1.0/rho)*px -mu/rho*dux + 1/rho*sux - mppf->chL/rho*lapPhi*phix;
      // mppf->rkU[id + 1*mppf->fieldOffset] = -nuy + (1.0/mppf->rho0 - 1.0/rho)*py -mu/rho*duy + 1/rho*suy - mppf->chL/rho*lapPhi*phiy;
       
      // dfloat t = time;
      // // Add forcing
      // dfloat gx = M_PI*cos(M_PI*x)*sin(M_PI*y)*cos(t) - sin(M_PI*x)*(cos(M_PI*x)*cos(M_PI*y)*sin(t) - 2.f)*(M_PI*cos(M_PI*x) + cos(M_PI*y)*cos(t) - M_PI*cos(M_PI*x)*cos(t)*cos(t)) - (M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*x)*sin(t)*sin(t))/100.f - 2.f*M_PI*M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(t)*((cos(M_PI*x)*cos(M_PI*y)*sin(t))/200.f - 3.f/200.f) + 2.f*M_PI*M_PI*M_PI*mppf->chL*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*x)*sin(t)*sin(t);
      // dfloat gy = sin(M_PI*y)*(cos(M_PI*x)*cos(M_PI*y)*sin(t) - 2.f)*(cos(M_PI*x)*cos(t) - M_PI*cos(M_PI*y) + M_PI*cos(M_PI*y)*cos(t)*cos(t)) + M_PI*cos(M_PI*y)*sin(M_PI*x)*cos(t) + (M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*y)*sin(t)*sin(t))/100.f + 2.f*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(t)*((cos(M_PI*x)*cos(M_PI*y)*sin(t))/200.f - 3.f/200.f) + 2.f*M_PI*M_PI*M_PI*mppf->chL*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*y)*sin(t)*sin(t);

      // mppf->rkU[id + 0*mppf->fieldOffset] += 1.0/rho*gx;
      // mppf->rkU[id + 1*mppf->fieldOffset] += 1.0/rho*gy;

      // mppf->rkU[id + 0*mppf->fieldOffset] *= mppf->dt;
      // mppf->rkU[id + 1*mppf->fieldOffset] *= mppf->dt;




    }
  }

  // o_rkU.copyFrom(mppf->rkU);
  
  // mppf->o_DU.copyFrom(mppf->DU); 
  // mppf->o_SU.copyFrom(mppf->SU); 
  // mppf->o_Rho.copyFrom(mppf->Rho); 
  // mppf->o_Mu.copyFrom(mppf->Mu); 
  // mppf->o_GP.copyFrom(mppf->GP); 
  
  // mppf->o_GPhi.copyFrom(mppf->GPhi); 
  // mppf->o_Psi.copyFrom(mppf->Psi); 
  // mppf->o_NU.copyFrom(mppf->NU); 

#endif




  // Upadte i.e. o_rkU = -NU + U^/dt + (1/rh0 - 1/rho^n+1)*GP^(*,n+1) -mu^(n+1)/rho^(n+1)*DU 
  //                                 + 1/rho^(n+1)*SU - lamdba*Psi*Ghi + 1/(rho^n+1)*f^(n+1) 
  mppf->explicitUpdateKernel(mesh->Nelements,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mppf->dt,
                           time,
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
  // o_rkU.copyTo(mppf->rkU);

  for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      dfloat rhsX =  cos(time)*(cos(M_PI*y)*sin(M_PI*x) + M_PI*cos(M_PI*x)*sin(M_PI*y));
      dfloat rhsY = -cos(time)*(cos(M_PI*x)*sin(M_PI*y) - M_PI*cos(M_PI*y)*sin(M_PI*x));

      mppf->rkU[id + 0*mppf->fieldOffset]  = mppf->dt*rhsX;
      mppf->rkU[id + 1*mppf->fieldOffset]  = mppf->dt*rhsY;
 
    }
  }

  // mppf->o_GSave.copyFrom(mppf->rkU); 
  o_rkU.copyFrom(mppf->rkU); 

#endif

  
  // rhsP = Div Uhat
  mppfDivergence(mppf, time, o_rkU, mppf->o_rhsP);




// Give exact NU
#if 0
  mppf->o_rhsP.copyTo(mppf->rhsP);
  for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      dfloat rhsp = -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(time);
 

     
      mppf->rhsP[id] -= mppf->dt*rhsp;
      // mppf->rhsP[id] = rhsp;
    }
  }

  // mppf->o_rhsP.copyFrom(mppf->rhsP); 
  mppf->o_Rho.copyFrom(mppf->rhsP); 

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
