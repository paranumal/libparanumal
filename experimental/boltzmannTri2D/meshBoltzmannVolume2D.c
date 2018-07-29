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

#include "mesh2D.h"

// function to compute collocation differentiation
// contributions to nodal DG rhs for boltzmann
void meshBoltzmannVolume2D(mesh2D *mesh){

  // for all elements
  for(int e=0;e<mesh->Nelements;++e){

    // prefetch geometric factors (constant on triangle)
    dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];

    // for all nodes in this element
    for(int n=0;n<mesh->Np;++n){
      // compute 'r' and 's' derivatives of (q_m) at node n
      dfloat dq1dr = 0, dq2dr = 0, dq3dr = 0, dq4dr = 0, dq5dr = 0, dq6dr = 0;
      dfloat dq1ds = 0, dq2ds = 0, dq3ds = 0, dq4ds = 0, dq5ds = 0, dq6ds = 0;
      
      for(int i=0;i<mesh->Np;++i){
	// load data at node i of element e (note Nfields==4)
	int id = mesh->Nfields*(e*mesh->Np + i);
	dfloat q1 = mesh->q[id+0];
	dfloat q2 = mesh->q[id+1];
	dfloat q3 = mesh->q[id+2];
	dfloat q4 = mesh->q[id+3];
	dfloat q5 = mesh->q[id+4];
	dfloat q6 = mesh->q[id+5];
	
	dfloat Drni = mesh->Dr[n*mesh->Np+i];
	dfloat Dsni = mesh->Ds[n*mesh->Np+i];

	// differentiate (u,v,p) with respect to 'r' and 's'
	dq1dr += Drni*q1;
	dq1ds += Dsni*q1;
	dq2dr += Drni*q2;
	dq2ds += Dsni*q2;
	dq3dr += Drni*q3;
	dq3ds += Dsni*q3;
	dq4dr += Drni*q4;
	dq4ds += Dsni*q4;
	dq5dr += Drni*q5;
	dq5ds += Dsni*q5;
	dq6dr += Drni*q6;
	dq6ds += Dsni*q6;
      }

      // chain rule
      dfloat dq1dx = drdx*dq1dr + dsdx*dq1ds;
      dfloat dq2dx = drdx*dq2dr + dsdx*dq2ds;
      dfloat dq3dx = drdx*dq3dr + dsdx*dq3ds;
      dfloat dq4dx = drdx*dq4dr + dsdx*dq4ds;
      dfloat dq5dx = drdx*dq5dr + dsdx*dq5ds;
      dfloat dq6dx = drdx*dq6dr + dsdx*dq6ds;

      dfloat dq1dy = drdy*dq1dr + dsdy*dq1ds;
      dfloat dq2dy = drdy*dq2dr + dsdy*dq2ds;
      dfloat dq3dy = drdy*dq3dr + dsdy*dq3ds;
      dfloat dq4dy = drdy*dq4dr + dsdy*dq4ds;
      dfloat dq5dy = drdy*dq5dr + dsdy*dq5ds;
      dfloat dq6dy = drdy*dq6dr + dsdy*dq6ds;
      
      // indices for writing the RHS terms
      int id = mesh->Nfields*(e*mesh->Np + n);

      // constant
      dfloat sqrtRT = mesh->sqrtRT, sqrt2 = sqrt(2.);
      dfloat invsqrt2 = 1./sqrt(2.);

      // transport operator
      dfloat rhsq1 = -sqrtRT*(dq2dx + dq3dy);
      dfloat rhsq2 = -sqrtRT*(dq1dx + sqrt2*dq5dx + dq4dy);	
      dfloat rhsq3 = -sqrtRT*(dq4dx + dq1dy + sqrt2*dq6dy);	
      dfloat rhsq4 = -sqrtRT*(dq3dx + dq2dy);		
      dfloat rhsq5 = -sqrtRT*sqrt2*dq2dx;			
      dfloat rhsq6 = -sqrtRT*sqrt2*dq3dy;                   

      // BGK relaxation approximation to the Boltzmann collision operator
      dfloat q1 = mesh->q[id+0];
      dfloat q2 = mesh->q[id+1];
      dfloat q3 = mesh->q[id+2];
      dfloat q4 = mesh->q[id+3];
      dfloat q5 = mesh->q[id+4];
      dfloat q6 = mesh->q[id+5];

      rhsq4 -= mesh->tauInv*(q4 - q2*q3/q1);
      rhsq5 -= mesh->tauInv*(q5 - invsqrt2*q2*q2/q1);
      rhsq6 -= mesh->tauInv*(q6 - invsqrt2*q3*q3/q1);
      
      // store boltzmann rhs contributions from collocation differentiation
      mesh->rhsq[id+0] = rhsq1;
      mesh->rhsq[id+1] = rhsq2;
      mesh->rhsq[id+2] = rhsq3;
      mesh->rhsq[id+3] = rhsq4;
      mesh->rhsq[id+4] = rhsq5;
      mesh->rhsq[id+5] = rhsq6;
    }
  }
}
