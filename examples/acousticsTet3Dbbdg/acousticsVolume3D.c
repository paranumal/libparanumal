#include "acoustics3D.h"



void acousticsVolume3Dbbdg(mesh3D *mesh, int lev){

  // for all elements
  for(int et=0;et<mesh->MRABNelements[lev];++et){
    int e = mesh->MRABelementIds[lev][et];

    // prefetch geometric factors (constant on tet)
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
      
      // compute 'r', 's', and 't' derivatives of (q_m) at node n
      int D0i1 = mesh->Nfields*(e*mesh->Np + mesh->D0ids[4*n+0]);
      int D1i1 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[4*n+0]);
      int D2i1 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[4*n+0]);
      int D3i1 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[4*n+0]);
      dfloat Dval1 = mesh->Dvals[4*n+0];
      
      int D0i2 = mesh->Nfields*(e*mesh->Np + mesh->D0ids[4*n+1]);
      int D1i2 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[4*n+1]);
      int D2i2 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[4*n+1]);
      int D3i2 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[4*n+1]);
      dfloat Dval2 = mesh->Dvals[4*n+1];

      int D0i3 = mesh->Nfields*(e*mesh->Np + mesh->D0ids[4*n+2]);
      int D1i3 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[4*n+2]);
      int D2i3 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[4*n+2]);
      int D3i3 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[4*n+2]);    
      dfloat Dval3 = mesh->Dvals[4*n+2];

      int D0i4 = mesh->Nfields*(e*mesh->Np + mesh->D0ids[4*n+3]);
      int D1i4 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[4*n+3]);
      int D2i4 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[4*n+3]);
      int D3i4 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[4*n+3]);    
      dfloat Dval4 = mesh->Dvals[4*n+3];

      dfloat dudr = .5f*(Dval1*(mesh->q[D1i1+0] - mesh->q[D0i1+0]) +
                         Dval2*(mesh->q[D1i2+0] - mesh->q[D0i2+0]) +
                         Dval3*(mesh->q[D1i3+0] - mesh->q[D0i3+0]) +
                         Dval4*(mesh->q[D1i4+0] - mesh->q[D0i4+0]));
      dfloat duds = .5f*(Dval1*(mesh->q[D2i1+0] - mesh->q[D0i1+0]) +
                         Dval2*(mesh->q[D2i2+0] - mesh->q[D0i2+0]) +
                         Dval3*(mesh->q[D2i3+0] - mesh->q[D0i3+0]) +
                         Dval4*(mesh->q[D2i4+0] - mesh->q[D0i4+0]));
      dfloat dudt = .5f*(Dval1*(mesh->q[D3i1+0] - mesh->q[D0i1+0]) +
                         Dval2*(mesh->q[D3i2+0] - mesh->q[D0i2+0]) +
                         Dval3*(mesh->q[D3i3+0] - mesh->q[D0i3+0]) +
                         Dval4*(mesh->q[D3i4+0] - mesh->q[D0i4+0]));

      dfloat dvdr = .5f*(Dval1*(mesh->q[D1i1+1] - mesh->q[D0i1+1]) +
                         Dval2*(mesh->q[D1i2+1] - mesh->q[D0i2+1]) +
                         Dval3*(mesh->q[D1i3+1] - mesh->q[D0i3+1]) +
                         Dval4*(mesh->q[D1i4+1] - mesh->q[D0i4+1]));
      dfloat dvds = .5f*(Dval1*(mesh->q[D2i1+1] - mesh->q[D0i1+1]) +
                         Dval2*(mesh->q[D2i2+1] - mesh->q[D0i2+1]) +
                         Dval3*(mesh->q[D2i3+1] - mesh->q[D0i3+1]) +
                         Dval4*(mesh->q[D2i4+1] - mesh->q[D0i4+1]));
      dfloat dvdt = .5f*(Dval1*(mesh->q[D3i1+1] - mesh->q[D0i1+1]) +
                         Dval2*(mesh->q[D3i2+1] - mesh->q[D0i2+1]) +
                         Dval3*(mesh->q[D3i3+1] - mesh->q[D0i3+1]) +
                         Dval4*(mesh->q[D3i4+1] - mesh->q[D0i4+1]));

      dfloat dwdr = .5f*(Dval1*(mesh->q[D1i1+2] - mesh->q[D0i1+2]) +
                         Dval2*(mesh->q[D1i2+2] - mesh->q[D0i2+2]) +
                         Dval3*(mesh->q[D1i3+2] - mesh->q[D0i3+2]) +
                         Dval4*(mesh->q[D1i4+2] - mesh->q[D0i4+2]));
      dfloat dwds = .5f*(Dval1*(mesh->q[D2i1+2] - mesh->q[D0i1+2]) +
                         Dval2*(mesh->q[D2i2+2] - mesh->q[D0i2+2]) +
                         Dval3*(mesh->q[D2i3+2] - mesh->q[D0i3+2]) +
                         Dval4*(mesh->q[D2i4+2] - mesh->q[D0i4+2]));
      dfloat dwdt = .5f*(Dval1*(mesh->q[D3i1+2] - mesh->q[D0i1+2]) +
                         Dval2*(mesh->q[D3i2+2] - mesh->q[D0i2+2]) +
                         Dval3*(mesh->q[D3i3+2] - mesh->q[D0i3+2]) +
                         Dval4*(mesh->q[D3i4+2] - mesh->q[D0i4+2]));

      dfloat dpdr = .5f*(Dval1*(mesh->q[D1i1+3] - mesh->q[D0i1+3]) +
                         Dval2*(mesh->q[D1i2+3] - mesh->q[D0i2+3]) +
                         Dval3*(mesh->q[D1i3+3] - mesh->q[D0i3+3]) +
                         Dval4*(mesh->q[D1i4+3] - mesh->q[D0i4+3]));
      dfloat dpds = .5f*(Dval1*(mesh->q[D2i1+3] - mesh->q[D0i1+3]) +
                         Dval2*(mesh->q[D2i2+3] - mesh->q[D0i2+3]) +
                         Dval3*(mesh->q[D2i3+3] - mesh->q[D0i3+3]) +
                         Dval4*(mesh->q[D2i4+3] - mesh->q[D0i4+3]));
      dfloat dpdt = .5f*(Dval1*(mesh->q[D3i1+3] - mesh->q[D0i1+3]) +
                         Dval2*(mesh->q[D3i2+3] - mesh->q[D0i2+3]) +
                         Dval3*(mesh->q[D3i3+3] - mesh->q[D0i3+3]) +
                         Dval4*(mesh->q[D3i4+3] - mesh->q[D0i4+3]));

      // chain rule
      dfloat dudx = drdx*dudr + dsdx*duds + dtdx*dudt;
      dfloat dvdy = drdy*dvdr + dsdy*dvds + dtdy*dvdt;
      dfloat dwdz = drdz*dwdr + dsdz*dwds + dtdz*dwdt;
      
      dfloat dpdx = drdx*dpdr + dsdx*dpds + dtdx*dpdt;
      dfloat dpdy = drdy*dpdr + dsdy*dpds + dtdy*dpdt;
      dfloat dpdz = drdz*dpdr + dsdz*dpds + dtdz*dpdt;
      
      // indices for writing the RHS terms
      int id = 3*mesh->Nfields*(e*mesh->Np + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      
      // store acoustics rhs contributions from collocation differentiation
      mesh->rhsq[id+0] = -dpdx;
      mesh->rhsq[id+1] = -dpdy;
      mesh->rhsq[id+2] = -dpdz;
      mesh->rhsq[id+3] = -dudx-dvdy-dwdz;
    }
  }
}