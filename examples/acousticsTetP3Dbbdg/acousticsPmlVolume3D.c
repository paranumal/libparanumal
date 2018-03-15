#include "acoustics3D.h"



void acousticsPmlVolume3Dbbdg(mesh3D *mesh, int lev){

  dfloat *cubu = (dfloat *) calloc(mesh->cubNpMax,sizeof(dfloat));
  dfloat *cubv = (dfloat *) calloc(mesh->cubNpMax,sizeof(dfloat));
  dfloat *cubw = (dfloat *) calloc(mesh->cubNpMax,sizeof(dfloat));

  dfloat *cubpx = (dfloat *) calloc(mesh->cubNpMax,sizeof(dfloat));
  dfloat *cubpy = (dfloat *) calloc(mesh->cubNpMax,sizeof(dfloat));
  dfloat *cubpz = (dfloat *) calloc(mesh->cubNpMax,sizeof(dfloat));

  // for all elements
  for(int et=0;et<mesh->MRABpmlNelements[lev];++et){
    int e = mesh->MRABpmlElementIds[lev][et];
    int pmlId = mesh->MRABpmlIds[lev][et];
    int N = mesh->N[e];

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
    for(int n=0;n<mesh->Np[N];++n){

      // compute 'r', 's', and 't' derivatives of (q_m) at node n
      int D0i1 = mesh->Nfields*(e*mesh->NpMax + mesh->D0ids[N][4*n+0]);
      int D1i1 = mesh->Nfields*(e*mesh->NpMax + mesh->D1ids[N][4*n+0]);
      int D2i1 = mesh->Nfields*(e*mesh->NpMax + mesh->D2ids[N][4*n+0]);
      int D3i1 = mesh->Nfields*(e*mesh->NpMax + mesh->D3ids[N][4*n+0]);
      dfloat Dval1 = mesh->Dvals[N][4*n+0];

      int D0i2 = mesh->Nfields*(e*mesh->NpMax + mesh->D0ids[N][4*n+1]);
      int D1i2 = mesh->Nfields*(e*mesh->NpMax + mesh->D1ids[N][4*n+1]);
      int D2i2 = mesh->Nfields*(e*mesh->NpMax + mesh->D2ids[N][4*n+1]);
      int D3i2 = mesh->Nfields*(e*mesh->NpMax + mesh->D3ids[N][4*n+1]);
      dfloat Dval2 = mesh->Dvals[N][4*n+1];

      int D0i3 = mesh->Nfields*(e*mesh->NpMax + mesh->D0ids[N][4*n+2]);
      int D1i3 = mesh->Nfields*(e*mesh->NpMax + mesh->D1ids[N][4*n+2]);
      int D2i3 = mesh->Nfields*(e*mesh->NpMax + mesh->D2ids[N][4*n+2]);
      int D3i3 = mesh->Nfields*(e*mesh->NpMax + mesh->D3ids[N][4*n+2]);
      dfloat Dval3 = mesh->Dvals[N][4*n+2];

      int D0i4 = mesh->Nfields*(e*mesh->NpMax + mesh->D0ids[N][4*n+3]);
      int D1i4 = mesh->Nfields*(e*mesh->NpMax + mesh->D1ids[N][4*n+3]);
      int D2i4 = mesh->Nfields*(e*mesh->NpMax + mesh->D2ids[N][4*n+3]);
      int D3i4 = mesh->Nfields*(e*mesh->NpMax + mesh->D3ids[N][4*n+3]);
      dfloat Dval4 = mesh->Dvals[N][4*n+3];

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
      int rhsid = 3*mesh->Nfields*(e*mesh->NpMax + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      int pmlrhsid = 3*mesh->pmlNfields*(pmlId*mesh->NpMax + n) + mesh->pmlNfields*mesh->MRABshiftIndex[lev];

      // store acoustics rhs contributions from collocation differentiation
      mesh->rhsq[rhsid+0] = -dpdx;
      mesh->rhsq[rhsid+1] = -dpdy;
      mesh->rhsq[rhsid+2] = -dpdz;

      mesh->pmlrhsq[pmlrhsid+0] = -dudx;
      mesh->pmlrhsq[pmlrhsid+1] = -dvdy;
      mesh->pmlrhsq[pmlrhsid+2] = -dwdz;
    }

    // Interpolate to cubature nodes
    for(int n=0;n<mesh->cubNp[N];++n){
      dfloat u = 0.f;
      dfloat v = 0.f;
      dfloat w = 0.f;
      dfloat px = 0.f;
      dfloat py = 0.f;
      dfloat pz = 0.f;
      for (int i=0;i<mesh->Np[N];++i){
        int base = mesh->Nfields*(e*mesh->NpMax + i);
        u += mesh->cubInterp[N][n*mesh->Np[N] + i] * mesh->q[base+0];
        v += mesh->cubInterp[N][n*mesh->Np[N] + i] * mesh->q[base+1];
        w += mesh->cubInterp[N][n*mesh->Np[N] + i] * mesh->q[base+2];

        int pmlBase = mesh->pmlNfields*(pmlId*mesh->NpMax+i);
        px += mesh->cubInterp[N][n*mesh->Np[N] + i] * mesh->pmlq[pmlBase+0];
        py += mesh->cubInterp[N][n*mesh->Np[N] + i] * mesh->pmlq[pmlBase+1];
        pz += mesh->cubInterp[N][n*mesh->Np[N] + i] * mesh->pmlq[pmlBase+2];
      }

      dfloat sigmax = mesh->pmlSigmaX[pmlId*mesh->cubNpMax + n];
      dfloat sigmay = mesh->pmlSigmaY[pmlId*mesh->cubNpMax + n];
      dfloat sigmaz = mesh->pmlSigmaZ[pmlId*mesh->cubNpMax + n];

      cubu[n] = -sigmax*u;
      cubv[n] = -sigmay*v;
      cubw[n] = -sigmaz*w;
      cubpx[n] = -sigmax*px;
      cubpy[n] = -sigmay*py;
      cubpz[n] = -sigmaz*pz;
    }

    //Project down and store
    for(int n=0;n<mesh->Np[N];++n){
      int rhsid = 3*mesh->Nfields*(e*mesh->NpMax + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      int pmlrhsid = 3*mesh->pmlNfields*(pmlId*mesh->NpMax + n) + mesh->pmlNfields*mesh->MRABshiftIndex[lev];

      dfloat rhsu = mesh->rhsq[rhsid+0];
      dfloat rhsv = mesh->rhsq[rhsid+1];
      dfloat rhsw = mesh->rhsq[rhsid+2];
      dfloat rhspx = mesh->pmlrhsq[pmlrhsid+0];
      dfloat rhspy = mesh->pmlrhsq[pmlrhsid+1];
      dfloat rhspz = mesh->pmlrhsq[pmlrhsid+2];
      for (int i=0;i<mesh->cubNp[N];++i){
        rhsu += mesh->cubProject[N][n*mesh->cubNp[N] + i] * cubu[i];
        rhsv += mesh->cubProject[N][n*mesh->cubNp[N] + i] * cubv[i];
        rhsw += mesh->cubProject[N][n*mesh->cubNp[N] + i] * cubw[i];
        rhspx += mesh->cubProject[N][n*mesh->cubNp[N] + i] * cubpx[i];
        rhspy += mesh->cubProject[N][n*mesh->cubNp[N] + i] * cubpy[i];
        rhspz += mesh->cubProject[N][n*mesh->cubNp[N] + i] * cubpz[i];
      }

      mesh->rhsq[rhsid+0] = rhsu;
      mesh->rhsq[rhsid+1] = rhsv;
      mesh->rhsq[rhsid+2] = rhsw;
      mesh->pmlrhsq[pmlrhsid+0] = rhspx;
      mesh->pmlrhsq[pmlrhsid+1] = rhspy;
      mesh->pmlrhsq[pmlrhsid+2] = rhspz;
    }
  }

  free(cubu); free(cubv); free(cubw);
  free(cubpx); free(cubpy); free(cubpz);
}