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

#include "bns.h"

void brownMinion(dfloat bmRho, dfloat bmDelta, dfloat sphereRadius,
		 dfloat x, dfloat y, dfloat z,
		 dfloat *u, dfloat *v, dfloat *w){

  dfloat Utangential = 0.25*(1+tanh(bmRho*(-z+0.5)))*(1+tanh(bmRho*(0.5+z)));

  dfloat uout, vout;

  if(x*x+y*y>1e-4) {
    uout =  -y*Utangential/(x*x+y*y);
    vout =   x*Utangential/(x*x+y*y);
  }
  else{
    uout = 0;
    vout = 0;
  }

  dfloat wout = bmDelta*sin(2*atan2(y,x))*(1-z*z);

  dfloat udotx = uout*x+vout*y+wout*z;
  *u = uout - udotx*x/(sphereRadius*sphereRadius);
  *v = vout - udotx*y/(sphereRadius*sphereRadius);
  *w = wout - udotx*z/(sphereRadius*sphereRadius);

}

void bnsBrownMinionQuad3D(bns_t *bns){

  mesh_t *mesh = bns->mesh;

  dfloat sR = mesh->sphereRadius;
  
  for(hlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){

      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];
      
      // Brown Minion shear layer roll up
      dfloat bmRho = 40;
      dfloat bmDelta  = 0.05;
      
      dfloat rho = 1;
  
      dfloat umod, vmod, wmod;
  
      brownMinion(bmRho, bmDelta, sR, x, y, z, &umod, &vmod, &wmod);
  
      dfloat delta = 1e-5;
  
      dfloat uP, uM, vP, vM, wP, wM;
  
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x+delta, y, z, &uP, &vP, &wP);
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x-delta, y, z, &uM, &vM, &wM);
  
      dfloat dudx = (uP-uM)/(2*delta);
      dfloat dvdx = (vP-vM)/(2*delta);
      dfloat dwdx = (wP-wM)/(2*delta);
  
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y+delta, z, &uP, &vP, &wP);
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y-delta, z, &uM, &vM, &wM);
  
      dfloat dudy = (uP-uM)/(2*delta);
      dfloat dvdy = (vP-vM)/(2*delta);
      dfloat dwdy = (wP-wM)/(2*delta);
  
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y, z+delta, &uP, &vP, &wP);
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y, z-delta, &uM, &vM, &wM);
  
      dfloat dudz = (uP-uM)/(2*delta);
      dfloat dvdz = (vP-vM)/(2*delta);
      dfloat dwdz = (wP-wM)/(2*delta);
  
      dfloat divu = dudx + dvdy + dwdz;
  
      dfloat sigma11 = bns->nu*(dudx+dudx - (2*divu/3));
      dfloat sigma12 = bns->nu*(dvdx+dudy);
      dfloat sigma13 = bns->nu*(dwdx+dudz);
      dfloat sigma22 = bns->nu*(dvdy+dvdy - (2*divu/3));
      dfloat sigma23 = bns->nu*(dwdy+dvdz);
      dfloat sigma33 = bns->nu*(dwdz+dwdz - (2*divu/3));

      dfloat q1bar = rho;
      dfloat q2bar = rho*umod/bns->sqrtRT;
      dfloat q3bar = rho*vmod/bns->sqrtRT;
      dfloat q4bar = rho*wmod/bns->sqrtRT;
      dfloat q5bar = (rho*umod*umod - sigma11)/(sqrt(2.)*bns->RT);
      dfloat q6bar = (rho*vmod*vmod - sigma22)/(sqrt(2.)*bns->RT);
      dfloat q7bar = (rho*wmod*wmod - sigma33)/(sqrt(2.)*bns->RT);
      dfloat q8bar  = (rho*umod*vmod - sigma12)/bns->RT;
      dfloat q9bar =  (rho*umod*wmod - sigma13)/bns->RT;
      dfloat q10bar = (rho*vmod*wmod - sigma23)/bns->RT;
  
      int base = n + e*mesh->Np*mesh->Nfields;
  
      bns->q[base+0*mesh->Np] = q1bar; // uniform density, zero flow
  
      bns->q[base+1*mesh->Np] = q2bar;
      bns->q[base+2*mesh->Np] = q3bar;
      bns->q[base+3*mesh->Np] = q4bar;
  
      bns->q[base+4*mesh->Np] = q5bar;
      bns->q[base+5*mesh->Np] = q6bar;
      bns->q[base+6*mesh->Np] = q7bar;
  
      bns->q[base+7*mesh->Np] = q8bar;
      bns->q[base+8*mesh->Np] = q9bar;
      bns->q[base+9*mesh->Np] = q10bar;
    }
  }
}
