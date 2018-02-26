#include "acoustics2D.h"

// function to compute collocation differentiation
// contributions to nodal DG rhs for acoustics
void acousticsVolume2D(mesh2D *mesh){

  // for all elements
  for(int e=0;e<mesh->Nelements;++e){

    // prefetch geometric factors (constant on triangle)
    dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];

    // for all nodes in this element
    for(int n=0;n<mesh->Np;++n){
      // compute 'r' and 's' derivatives of (u,v,p) at node n
      dfloat dudr = 0, duds = 0;
      dfloat dvdr = 0, dvds = 0;
      dfloat dpdr = 0, dpds = 0;
      
      for(int i=0;i<mesh->Np;++i){
	// load data at node i of element e
	int id = mesh->Nfields*(e*mesh->Np + i);
	dfloat u = mesh->q[id+0];
	dfloat v = mesh->q[id+1];
	dfloat p = mesh->q[id+2];

	dfloat Drni = mesh->Dr[n*mesh->Np+i];
	dfloat Dsni = mesh->Ds[n*mesh->Np+i];

	// differentiate (u,v,p) with respect to 'r' and 's'
	dudr += Drni*u;
	duds += Dsni*u;
	dvdr += Drni*v;
	dvds += Dsni*v;
	dpdr += Drni*p;
	dpds += Dsni*p;
      }

      // chain rule
      dfloat dudx = drdx*dudr + dsdx*duds;
      dfloat dvdy = drdy*dvdr + dsdy*dvds;
      dfloat dpdx = drdx*dpdr + dsdx*dpds;
      dfloat dpdy = drdy*dpdr + dsdy*dpds;
      
      // indices for writing the RHS terms
      int id = mesh->Nfields*(e*mesh->Np + n);

      // store acoustics rhs contributions from collocation differentiation
      mesh->rhsq[id+0] = -dpdx;
      mesh->rhsq[id+1] = -dpdy;
      mesh->rhsq[id+2] = -dudx-dvdy;
    }
  }
}
