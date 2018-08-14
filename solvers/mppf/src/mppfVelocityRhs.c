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

void mppfVelocityRhs(mppf_t *mppf, dfloat time, occa::memory o_rhsU, occa::memory o_rhsV, occa::memory o_rhsW){
  
  mesh_t *mesh = mppf->mesh; 

// Give exact DU and SU are fine
#if 0

  mppf->o_GP.copyTo(mppf->GP);


  for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      dfloat u     =  cos(M_PI*y)*sin(M_PI*x)*sin(time);
      dfloat v     = -sin(M_PI*y)*cos(M_PI*x)*sin(time);

      dfloat dux =  2.0*M_PI*M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(time);
      dfloat duy = -2.0*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(time);

      dfloat rho = 2.0 - cos(M_PI*x)*cos(M_PI*y)*sin(time);
      dfloat mu  = 3.0/200.0 - (cos(M_PI*x)*cos(M_PI*y)*sin(time))/200.0;

      dfloat sux =  (M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*x)*sin(time)*sin(time))/100.0;
      dfloat suy = -(M_PI*M_PI*cos(M_PI*x)*cos(M_PI*x)*cos(M_PI*y)*sin(M_PI*y)*sin(time)*sin(time))/100.0;

      dfloat px = M_PI*cos(M_PI*x)*sin(M_PI*y)*cos(time);
      dfloat py = M_PI*cos(M_PI*y)*sin(M_PI*x)*cos(time);

      dfloat lapPhi = -2.0*M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(time);
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



      mppf->rkU[id + 0*mppf->fieldOffset] = u/mppf->dt -nux + (1.0/mppf->rho0 - 1.0/rho)*px -mu/rho*dux + 1/rho*sux - mppf->chL/rho*lapPhi*phix;
      mppf->rkU[id + 1*mppf->fieldOffset] = v/mppf->dt -nuy + (1.0/mppf->rho0 - 1.0/rho)*py -mu/rho*duy + 1/rho*suy - mppf->chL/rho*lapPhi*phiy;
       
       dfloat t = time;
      // // Add forcing
      dfloat gx = M_PI*cos(M_PI*x)*sin(M_PI*y)*cos(t) - sin(M_PI*x)*(cos(M_PI*x)*cos(M_PI*y)*sin(t) - 2.f)*(M_PI*cos(M_PI*x) + cos(M_PI*y)*cos(t) - M_PI*cos(M_PI*x)*cos(t)*cos(t)) - (M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*x)*sin(t)*sin(t))/100.f - 2.f*M_PI*M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(t)*((cos(M_PI*x)*cos(M_PI*y)*sin(t))/200.f - 3.f/200.f) + 2.f*M_PI*M_PI*M_PI*mppf->chL*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*x)*sin(t)*sin(t);
      dfloat gy = sin(M_PI*y)*(cos(M_PI*x)*cos(M_PI*y)*sin(t) - 2.f)*(cos(M_PI*x)*cos(t) - M_PI*cos(M_PI*y) + M_PI*cos(M_PI*y)*cos(t)*cos(t)) + M_PI*cos(M_PI*y)*sin(M_PI*x)*cos(t) + (M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*y)*sin(t)*sin(t))/100.f + 2.f*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(t)*((cos(M_PI*x)*cos(M_PI*y)*sin(t))/200.f - 3.f/200.f) + 2.f*M_PI*M_PI*M_PI*mppf->chL*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*y)*sin(M_PI*y)*sin(t)*sin(t);

      mppf->rkU[id + 0*mppf->fieldOffset] += 1.0/rho*gx;
      mppf->rkU[id + 1*mppf->fieldOffset] += 1.0/rho*gy;

      mppf->rkU[id + 0*mppf->fieldOffset] *= mppf->dt;
      mppf->rkU[id + 1*mppf->fieldOffset] *= mppf->dt;




    }
  }

  mppf->o_Uhat.copyFrom(mppf->rkU);
  
  mppf->o_DU.copyFrom(mppf->DU); 
  // mppf->o_SU.copyFrom(mppf->SU); 
  // mppf->o_Rho.copyFrom(mppf->Rho); 
  // mppf->o_Mu.copyFrom(mppf->Mu); 
  mppf->o_GP.copyFrom(mppf->GP); 
  
  // mppf->o_GPhi.copyFrom(mppf->GPhi); 
  // mppf->o_Psi.copyFrom(mppf->Psi); 
  // mppf->o_NU.copyFrom(mppf->NU); 

#endif



     // rhsU^s = MM*(\sum^s b_i U^n-i - \sum^s-1 a_i N(U^n-i) + \sum^s-1 c_i GP^n-i)/nu dt
    mppf->velocityRhsKernel(mesh->Nelements,
                           mesh->o_vgeo,
                           mesh->o_MM,
                           mppf->idt,
                           mppf->fieldOffset,
                           mppf->o_Uhat,
                           mppf->o_DU,
                           mppf->o_GP,
                           o_rhsU,
                           o_rhsV,
                           o_rhsW);
}
