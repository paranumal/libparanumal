#include "acoustics3D.h"



void acousticsPmlVolume3Dbbdg(mesh3D *mesh, iint lev){

  dfloat *cubu = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubv = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubw = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));

  dfloat *cubpx = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubpy = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubpz = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));

  // for all elements
  for(iint et=0;et<mesh->MRABpmlNelements[lev];++et){
    iint e = mesh->MRABpmlElementIds[lev][et];
    iint pmlId = mesh->MRABpmlIds[lev][et];

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
    for(iint n=0;n<mesh->Np;++n){

      // compute 'r', 's', and 't' derivatives of (q_m) at node n
      iint D0i1 = mesh->Nfields*(e*mesh->Np + mesh->D0ids[4*n+0]);
      iint D1i1 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[4*n+0]);
      iint D2i1 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[4*n+0]);
      iint D3i1 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[4*n+0]);
      dfloat Dval1 = mesh->Dvals[4*n+0];

      iint D0i2 = mesh->Nfields*(e*mesh->Np + mesh->D0ids[4*n+1]);
      iint D1i2 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[4*n+1]);
      iint D2i2 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[4*n+1]);
      iint D3i2 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[4*n+1]);
      dfloat Dval2 = mesh->Dvals[4*n+1];

      iint D0i3 = mesh->Nfields*(e*mesh->Np + mesh->D0ids[4*n+2]);
      iint D1i3 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[4*n+2]);
      iint D2i3 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[4*n+2]);
      iint D3i3 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[4*n+2]);
      dfloat Dval3 = mesh->Dvals[4*n+2];

      iint D0i4 = mesh->Nfields*(e*mesh->Np + mesh->D0ids[4*n+3]);
      iint D1i4 = mesh->Nfields*(e*mesh->Np + mesh->D1ids[4*n+3]);
      iint D2i4 = mesh->Nfields*(e*mesh->Np + mesh->D2ids[4*n+3]);
      iint D3i4 = mesh->Nfields*(e*mesh->Np + mesh->D3ids[4*n+3]);
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
      iint rhsid = 3*mesh->Nfields*(e*mesh->Np + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      iint pmlrhsid = 3*mesh->pmlNfields*(pmlId*mesh->Np + n) + mesh->pmlNfields*mesh->MRABshiftIndex[lev];

      // store acoustics rhs contributions from collocation differentiation
      mesh->rhsq[rhsid+0] = -dpdx;
      mesh->rhsq[rhsid+1] = -dpdy;
      mesh->rhsq[rhsid+2] = -dpdz;

      mesh->pmlrhsq[pmlrhsid+0] = -dudx;
      mesh->pmlrhsq[pmlrhsid+1] = -dvdy;
      mesh->pmlrhsq[pmlrhsid+2] = -dwdz;
    }

    // Interpolate to cubature nodes
    for(iint n=0;n<mesh->cubNp;++n){
      dfloat u = 0.f;
      dfloat v = 0.f;
      dfloat w = 0.f;
      dfloat px = 0.f;
      dfloat py = 0.f;
      dfloat pz = 0.f;
      for (iint i=0;i<mesh->Np;++i){
        iint base = mesh->Nfields*(e*mesh->Np + i);
        u += mesh->cubInterp[n*mesh->Np + i] * mesh->q[base+0];
        v += mesh->cubInterp[n*mesh->Np + i] * mesh->q[base+1];
        w += mesh->cubInterp[n*mesh->Np + i] * mesh->q[base+2];

        iint pmlBase = mesh->pmlNfields*(pmlId*mesh->Np+i);
        px += mesh->cubInterp[n*mesh->Np + i] * mesh->pmlq[pmlBase+0];
        py += mesh->cubInterp[n*mesh->Np + i] * mesh->pmlq[pmlBase+1];
        pz += mesh->cubInterp[n*mesh->Np + i] * mesh->pmlq[pmlBase+2];
      }

      dfloat sigmax = mesh->pmlSigmaX[pmlId*mesh->cubNp + n];
      dfloat sigmay = mesh->pmlSigmaY[pmlId*mesh->cubNp + n];
      dfloat sigmaz = mesh->pmlSigmaZ[pmlId*mesh->cubNp + n];

      cubu[n] = -sigmax*u;
      cubv[n] = -sigmay*v;
      cubw[n] = -sigmaz*w;
      cubpx[n] = -sigmax*px;
      cubpy[n] = -sigmay*py;
      cubpz[n] = -sigmaz*pz;
    }

    //Project down and store
    for(iint n=0;n<mesh->Np;++n){
      iint rhsid = 3*mesh->Nfields*(e*mesh->Np + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      iint pmlrhsid = 3*mesh->pmlNfields*(pmlId*mesh->Np + n) + mesh->pmlNfields*mesh->MRABshiftIndex[lev];

      dfloat rhsu = mesh->rhsq[rhsid+0];
      dfloat rhsv = mesh->rhsq[rhsid+1];
      dfloat rhsw = mesh->rhsq[rhsid+2];
      dfloat rhspx = mesh->pmlrhsq[pmlrhsid+0];
      dfloat rhspy = mesh->pmlrhsq[pmlrhsid+1];
      dfloat rhspz = mesh->pmlrhsq[pmlrhsid+2];
      for (iint i=0;i<mesh->cubNp;++i){
        rhsu += mesh->cubProject[n*mesh->cubNp + i] * cubu[i];
        rhsv += mesh->cubProject[n*mesh->cubNp + i] * cubv[i];
        rhsw += mesh->cubProject[n*mesh->cubNp + i] * cubw[i];
        rhspx += mesh->cubProject[n*mesh->cubNp + i] * cubpx[i];
        rhspy += mesh->cubProject[n*mesh->cubNp + i] * cubpy[i];
        rhspz += mesh->cubProject[n*mesh->cubNp + i] * cubpz[i];
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