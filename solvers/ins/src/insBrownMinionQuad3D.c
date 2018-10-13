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

  //  dfloat wout = bmDelta*sin(2*atan2(y,x))*(1-z*z);
  dfloat wout = bmDelta*sin(8*atan2(y,x))*(1-z*z);

  dfloat udotx = uout*x+vout*y+wout*z;
  *u = uout - udotx*x/(sphereRadius*sphereRadius);
  *v = vout - udotx*y/(sphereRadius*sphereRadius);
  *w = wout - udotx*z/(sphereRadius*sphereRadius);

}

void insBrownMinionQuad3D(ins_t *ins){

  mesh_t *mesh = ins->mesh;

  dfloat sR = mesh->sphereRadius;
  
  for(hlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){

      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];
      
      // Brown Minion shear layer roll up
      dfloat bmRho    = 10.0;
      dfloat bmDelta  = 0.05;
      
      dfloat pr = 0.0;
  
      dfloat umod, vmod, wmod;
  
      brownMinion(bmRho, bmDelta, sR, x, y, z, &umod, &vmod, &wmod);
  
      int base = n + e*mesh->Np;
  
      ins->U[base+0*ins->fieldOffset] = umod; 
      ins->U[base+1*ins->fieldOffset] = vmod;
      ins->U[base+2*ins->fieldOffset] = wmod;
      
      ins->P[base+0*ins->fieldOffset] = pr;
    }
  }
}
