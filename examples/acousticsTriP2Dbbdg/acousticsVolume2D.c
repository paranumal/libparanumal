#include "acoustics2D.h"

void acousticsVolume2Dbbdg(mesh2D *mesh, int lev){

  // loop over elements
  for(int et=0;et<mesh->MRABNelements[lev];++et){
    int e = mesh->MRABelementIds[lev][et];
    
    // prefetch geometric factors (constant on triangle)
    dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];

    int N = mesh->N[e];

    for(int n=0;n<mesh->Np[N];++n){     // for all nodes in this element

      // compute 'r' and 's' derivatives of (q_m) at node n
      int D1i1 = mesh->Nfields*(e*mesh->NpMax + mesh->D1ids[N][3*n]);
      int D2i1 = mesh->Nfields*(e*mesh->NpMax + mesh->D2ids[N][3*n]);
      int D3i1 = mesh->Nfields*(e*mesh->NpMax + mesh->D3ids[N][3*n]);
      dfloat Dval1 = mesh->Dvals[N][3*n];
      
      int D1i2 = mesh->Nfields*(e*mesh->NpMax + mesh->D1ids[N][3*n+1]);
      int D2i2 = mesh->Nfields*(e*mesh->NpMax + mesh->D2ids[N][3*n+1]);
      int D3i2 = mesh->Nfields*(e*mesh->NpMax + mesh->D3ids[N][3*n+1]);
      dfloat Dval2 = mesh->Dvals[N][3*n+1];

      int D1i3 = mesh->Nfields*(e*mesh->NpMax + mesh->D1ids[N][3*n+2]);
      int D2i3 = mesh->Nfields*(e*mesh->NpMax + mesh->D2ids[N][3*n+2]);
      int D3i3 = mesh->Nfields*(e*mesh->NpMax + mesh->D3ids[N][3*n+2]);    
      dfloat Dval3 = mesh->Dvals[N][3*n+2];

      dfloat dudr = .5f*(Dval1*(mesh->q[D2i1+0] - mesh->q[D1i1+0]) +
                         Dval2*(mesh->q[D2i2+0] - mesh->q[D1i2+0]) +
                         Dval3*(mesh->q[D2i3+0] - mesh->q[D1i3+0]));
      dfloat duds = .5f*(Dval1*(mesh->q[D3i1+0] - mesh->q[D1i1+0]) +
                         Dval2*(mesh->q[D3i2+0] - mesh->q[D1i2+0]) +
                         Dval3*(mesh->q[D3i3+0] - mesh->q[D1i3+0]));

      dfloat dvdr = .5f*(Dval1*(mesh->q[D2i1+1] - mesh->q[D1i1+1]) +
                         Dval2*(mesh->q[D2i2+1] - mesh->q[D1i2+1]) +
                         Dval3*(mesh->q[D2i3+1] - mesh->q[D1i3+1]));
      dfloat dvds = .5f*(Dval1*(mesh->q[D3i1+1] - mesh->q[D1i1+1]) +
                         Dval2*(mesh->q[D3i2+1] - mesh->q[D1i2+1]) +
                         Dval3*(mesh->q[D3i3+1] - mesh->q[D1i3+1]));
      
      dfloat dpdr = .5f*(Dval1*(mesh->q[D2i1+2] - mesh->q[D1i1+2]) +
                         Dval2*(mesh->q[D2i2+2] - mesh->q[D1i2+2]) +
                         Dval3*(mesh->q[D2i3+2] - mesh->q[D1i3+2]));
      dfloat dpds = .5f*(Dval1*(mesh->q[D3i1+2] - mesh->q[D1i1+2]) +
                         Dval2*(mesh->q[D3i2+2] - mesh->q[D1i2+2]) +
                         Dval3*(mesh->q[D3i3+2] - mesh->q[D1i3+2]));

      // chain rule
      dfloat dudx = drdx*dudr + dsdx*duds;
      dfloat dvdy = drdy*dvdr + dsdy*dvds;
      dfloat dpdx = drdx*dpdr + dsdx*dpds;
      dfloat dpdy = drdy*dpdr + dsdy*dpds;
      
      // indices for writing the RHS terms
      int id = 3*mesh->Nfields*(e*mesh->NpMax + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];

      // store acoustics rhs contributions from collocation differentiation
      mesh->rhsq[id+0] = -dpdx;
      mesh->rhsq[id+1] = -dpdy;
      mesh->rhsq[id+2] = -dudx-dvdy;
    }
  }
}