#include "cnsQuad2D.h"

// function to compute collocation differentiation
// contributions to nodal DG rhs for cns
void cnsVolumeQuad2D(mesh2D *mesh){

  // for all elements
  for(int e=0;e<mesh->Nelements;++e){

    // for all nodes in this element
    for(int m=0;m<mesh->Nq;++m){
      for(int n=0;n<mesh->Nq;++n){

	int id = n + m*mesh->Nq;
	
	// prefetch geometric factors (constant on triangle)
	dfloat drdx = mesh->vgeo[id+mesh->Np*(e*mesh->Nvgeo + RXID)];
	dfloat drdy = mesh->vgeo[id+mesh->Np*(e*mesh->Nvgeo + RYID)];
	dfloat dsdx = mesh->vgeo[id+mesh->Np*(e*mesh->Nvgeo + SXID)];
	dfloat dsdy = mesh->vgeo[id+mesh->Np*(e*mesh->Nvgeo + SYID)];
      
	// compute 'r' and 's' derivatives of (u,v,p) at node n
	dfloat dudr = 0, duds = 0;
	dfloat dvdr = 0, dvds = 0;
	dfloat dpdr = 0, dpds = 0;

	// 'r' derivatives
	for(int i=0;i<mesh->Nq;++i){

	  int id = mesh->Nfields*(e*mesh->Np + i + m*mesh->Nq);
	  dfloat u = mesh->q[id+0];
	  dfloat v = mesh->q[id+1];
	  dfloat p = mesh->q[id+2];
	  
	  dfloat Dni = mesh->D[n*mesh->Nq+i];

	  dudr += Dni*u;
	  dvdr += Dni*v;
	  dpdr += Dni*p;
	}

	// 's' derivatives
	for(int i=0;i<mesh->Nq;++i){

	  int id = mesh->Nfields*(e*mesh->Np + n + i*mesh->Nq);
	  dfloat u = mesh->q[id+0];
	  dfloat v = mesh->q[id+1];
	  dfloat p = mesh->q[id+2];
	  
	  dfloat Dmi = mesh->D[m*mesh->Nq+i];

	  duds += Dmi*u;
	  dvds += Dmi*v;
	  dpds += Dmi*p;
	}

	// chain rule
	dfloat dudx = drdx*dudr + dsdx*duds;
	dfloat dvdy = drdy*dvdr + dsdy*dvds;
	dfloat dpdx = drdx*dpdr + dsdx*dpds;
	dfloat dpdy = drdy*dpdr + dsdy*dpds;
	
	// indices for writing the RHS terms
	id = mesh->Nfields*(e*mesh->Np + n + m*mesh->Nq);
	
	// store cns rhs contributions from collocation differentiation
	mesh->rhsq[id+0] = -dpdx;
	mesh->rhsq[id+1] = -dpdy;
	mesh->rhsq[id+2] = -dudx-dvdy;
      }
    }
  }
}
