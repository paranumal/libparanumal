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

#include "acoustics3D.h"

// function to compute collocation differentiation
// contributions to nodal DG rhs for acoustics
void acousticsVolume3D(mesh3D *mesh){

  // for all elements
  for(int e=0;e<mesh->Nelements;++e){

    // prefetch geometric factors (constant on triangle)
    dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat drdz = mesh->vgeo[e*mesh->Nvgeo + RZID];
    dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];
    dfloat dsdz = mesh->vgeo[e*mesh->Nvgeo + SZID];
    dfloat dtdx = mesh->vgeo[e*mesh->Nvgeo + TXID];
    dfloat dtdy = mesh->vgeo[e*mesh->Nvgeo + TYID];
    dfloat dtdz = mesh->vgeo[e*mesh->Nvgeo + TZID];

    // for all nodes in this element
    for(int n=0;n<mesh->Np;++n){
      // compute 'r' and 's' derivatives of (u,v,p) at node n
      dfloat dudr = 0, duds = 0, dudt = 0;
      dfloat dvdr = 0, dvds = 0, dvdt = 0;
      dfloat dwdr = 0, dwds = 0, dwdt = 0;
      dfloat dpdr = 0, dpds = 0, dpdt = 0;
      
      for(int i=0;i<mesh->Np;++i){
	// load data at node i of element e
	int id = mesh->Nfields*(e*mesh->Np + i);
	dfloat u = mesh->q[id+0];
	dfloat v = mesh->q[id+1];
	dfloat w = mesh->q[id+2];
	dfloat p = mesh->q[id+3];

	dfloat Drni = mesh->Dr[n*mesh->Np+i];
	dfloat Dsni = mesh->Ds[n*mesh->Np+i];
	dfloat Dtni = mesh->Dt[n*mesh->Np+i];

	// differentiate (u,v,w,p) with respect to 'r,s,t'
	dudr += Drni*u;	duds += Dsni*u; dudt += Dtni*u;
	dvdr += Drni*v;	dvds += Dsni*v; dvdt += Dtni*v;
	dwdr += Drni*w;	dwds += Dsni*w; dwdt += Dtni*w;
	dpdr += Drni*p;	dpds += Dsni*p; dpdt += Dtni*p;
      }

      // chain rule
      dfloat dudx = drdx*dudr + dsdx*duds + dtdx*dudt;
      dfloat dvdy = drdy*dvdr + dsdy*dvds + dtdy*dvdt;
      dfloat dwdz = drdz*dwdr + dsdz*dwds + dtdz*dwdt;
      
      dfloat dpdx = drdx*dpdr + dsdx*dpds + dtdx*dpdt;
      dfloat dpdy = drdy*dpdr + dsdy*dpds + dtdy*dpdt;
      dfloat dpdz = drdz*dpdr + dsdz*dpds + dtdz*dpdt;
      
      // indices for writing the RHS terms
      int id = mesh->Nfields*(e*mesh->Np + n);

      // store acoustics rhs contributions from collocation differentiation
      mesh->rhsq[id+0] = -dpdx;
      mesh->rhsq[id+1] = -dpdy;
      mesh->rhsq[id+2] = -dpdz;
      mesh->rhsq[id+3] = -dudx-dvdy-dwdz;
    }
  }
}
