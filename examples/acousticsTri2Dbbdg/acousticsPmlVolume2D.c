#include "mesh2D.h"


void acousticsPmlVolume2Dbbdg(mesh2D *mesh, iint lev){

  // loop over elements
  for(iint et=0;et<mesh->MRABpmlNelements[lev];++et){
    iint e = mesh->MRABpmlElementIds[lev][et];
    iint pmlId = mesh->MRABpmlIds[lev][et];

    // prefetch geometric factors (constant on triangle)
    dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];

    for(iint n=0;n<mesh->Np;++n){     // for all nodes in this element

      dfloat sigmax = mesh->pmlSigmaX[pmlId*mesh->Np + n];
      dfloat sigmay = mesh->pmlSigmaY[pmlId*mesh->Np + n];

      iint base = mesh->Nfields*(e*mesh->Np+n);
      dfloat u = mesh->q[base+0];
      dfloat v = mesh->q[base+1];

      iint pmlBase = mesh->pmlNfields*(pmlId*mesh->Np+n);
      dfloat px = mesh->pmlq[pmlBase+0];
      dfloat py = mesh->pmlq[pmlBase+1];

      // compute 'r' and 's' derivatives of (q_m) at node n
      iint D1i1 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[3*n]);
      iint D2i1 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[3*n]);
      iint D3i1 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[3*n]);
      dfloat Dval1 = mesh->Dvals[3*n];
      
      iint D1i2 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[3*n+1]);
      iint D2i2 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[3*n+1]);
      iint D3i2 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[3*n+1]);
      dfloat Dval2 = mesh->Dvals[3*n+1];

      iint D1i3 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[3*n+2]);
      iint D2i3 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[3*n+2]);
      iint D3i3 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[3*n+2]);    
      dfloat Dval3 = mesh->Dvals[3*n+2];

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
      
      // chain rule
      dfloat dudx = drdx*dudr + dsdx*duds;
      dfloat dvdy = drdy*dvdr + dsdy*dvds;
      
      // indices for writing the RHS terms
      iint rhsId = 3*mesh->Nfields*(e*mesh->Np + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];

      // store acoustics pml rhs contributions 
      mesh->rhsq[rhsId+0] += -sigmax*u;
      mesh->rhsq[rhsId+1] += -sigmay*v;
      mesh->rhsq[rhsId+2] += -sigmax*px-sigmay*py;

      iint pmlrhsId = 3*mesh->pmlNfields*(pmlId*mesh->Np + n) + mesh->pmlNfields*mesh->MRABshiftIndex[lev];
      //store pml rhs
      mesh->pmlrhsq[pmlrhsId+0] = -dudx-sigmax*px;
      mesh->pmlrhsq[pmlrhsId+1] = -dvdy-sigmay*py;
    }
  }
}
