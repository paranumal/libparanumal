#include "acoustics2D.h"

// function to compute collocation differentiation
// contributions to nodal DG rhs for acoustics
void acousticsPml2D(mesh2D *mesh, int lev){

  dfloat *p = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *u = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *v = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));

  dfloat *ptilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *utilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *vtilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));

  dfloat *cubrhsp = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubrhsu = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubrhsv = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));

  dfloat *cubrhsptilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubrhsutilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubrhsvtilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));

  // for all elements
  for(iint m=0;m<mesh->MRABpmlNelements[lev];++m){
    iint e = mesh->MRABpmlElementIds[lev][m];
    iint pmlId = mesh->MRABpmlIds[lev][m];
    
    // Interpolate to cubature nodes
    for(iint n=0;n<mesh->cubNp;++n){
      p[n] = 0.f;
      u[n] = 0.f;
      v[n] = 0.f;
      ptilde[n] = 0.f;
      utilde[n] = 0.f;
      vtilde[n] = 0.f;
      for (iint i=0;i<mesh->Np;++i){
        iint base = mesh->Nfields*(e*mesh->Np + i);
        u[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->q[base+0];
        v[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->q[base+1];
        p[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->q[base+2];

        iint pmlBase = mesh->pmlNfields*(pmlId*mesh->Np+i);
        utilde[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->pmlq[pmlBase+0];
        vtilde[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->pmlq[pmlBase+1];
        ptilde[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->pmlq[pmlBase+2];        
      }    

      dfloat sigmax = mesh->pmlSigmaX[pmlId*mesh->cubNp + n];
      dfloat sigmay = mesh->pmlSigmaY[pmlId*mesh->cubNp + n];

      cubrhsu[n] = -(sigmax-sigmay)*u[n] - sigmay*(sigmay-sigmax)*utilde[n]; // uhat
      cubrhsv[n] = -(sigmay-sigmax)*v[n] - sigmax*(sigmax-sigmay)*vtilde[n]; // vhat
      cubrhsp[n] = -(sigmax+sigmay)*p[n] - sigmax*sigmay*ptilde[n]; // p

      // update for u~,v~, p~
      cubrhsutilde[n] = u[n]-sigmay*utilde[n]; // du~/dt = -sigmay*u~  + uhat
      cubrhsvtilde[n] = v[n]-sigmax*vtilde[n]; // dv~/dt = -sigmax*v~  + vhat
      cubrhsptilde[n] = p[n];                  // dp~/dt = p
    }

    //Project down and store
    for(iint n=0;n<mesh->Np;++n){
      dfloat rhsp = 0.f;
      dfloat rhsu = 0.f;
      dfloat rhsv = 0.f;
      dfloat rhsptilde = 0.f;
      dfloat rhsutilde = 0.f;
      dfloat rhsvtilde = 0.f;
      for (iint i=0;i<mesh->cubNp;++i){
        rhsp += mesh->cubProject[n*mesh->cubNp + i] * cubrhsp[i];
        rhsu += mesh->cubProject[n*mesh->cubNp + i] * cubrhsu[i];
        rhsv += mesh->cubProject[n*mesh->cubNp + i] * cubrhsv[i];
        rhsptilde += mesh->cubProject[n*mesh->cubNp + i] * cubrhsptilde[i];
        rhsutilde += mesh->cubProject[n*mesh->cubNp + i] * cubrhsutilde[i];
        rhsvtilde += mesh->cubProject[n*mesh->cubNp + i] * cubrhsvtilde[i];
      }

      // store acoustics rhs contributions from PML
      iint rhsBase = 3*mesh->Nfields*(e*mesh->Np + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      mesh->rhsq[rhsBase+0] += rhsu;
      mesh->rhsq[rhsBase+1] += rhsv;
      mesh->rhsq[rhsBase+2] += rhsp;

      // store PML rhs
      iint pmlrhsBase = 3*mesh->pmlNfields*(n+pmlId*mesh->Np)  + mesh->pmlNfields*mesh->MRABshiftIndex[lev];
      mesh->pmlrhsq[pmlrhsBase+0] = rhsutilde;
      mesh->pmlrhsq[pmlrhsBase+1] = rhsvtilde;
      mesh->pmlrhsq[pmlrhsBase+2] = rhsptilde;
    }
  }

  free(p); free(u); free(v);
  free(cubrhsp); free(cubrhsu); free(cubrhsv);
  free(ptilde); free(utilde); free(vtilde);
  free(cubrhsptilde); free(cubrhsutilde); free(cubrhsvtilde);
}


// function to compute collocation differentiation
// contributions to nodal DG rhs for acoustics
void acousticsPml2D_wadg(mesh2D *mesh, int lev){

  dfloat *p = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *u = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *v = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));

  dfloat *ptilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *utilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *vtilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));

  dfloat *cubrhsp = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubrhsu = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubrhsv = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));

  dfloat *cubrhsptilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubrhsutilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  dfloat *cubrhsvtilde = (dfloat *) calloc(mesh->cubNp,sizeof(dfloat));
  
  // for all elements
  for(iint m=0;m<mesh->MRABpmlNelements[lev];++m){
    iint e = mesh->MRABpmlElementIds[lev][m];
    iint pmlId = mesh->MRABpmlIds[lev][m];
    
    // Interpolate to cubature nodes
    for(iint n=0;n<mesh->cubNp;++n){
      p[n] = 0.f;
      u[n] = 0.f;
      v[n] = 0.f;
      ptilde[n] = 0.f;
      utilde[n] = 0.f;
      vtilde[n] = 0.f;
      for (iint i=0;i<mesh->Np;++i){
        iint base = mesh->Nfields*(e*mesh->Np + i);
        u[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->q[base+0];
        v[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->q[base+1];
        p[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->q[base+2];

        iint pmlBase = mesh->pmlNfields*(pmlId*mesh->Np+i);
        utilde[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->pmlq[pmlBase+0];
        vtilde[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->pmlq[pmlBase+1];
        ptilde[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->pmlq[pmlBase+2];        
      }    

      dfloat sigmax = mesh->pmlSigmaX[pmlId*mesh->cubNp + n];
      dfloat sigmay = mesh->pmlSigmaY[pmlId*mesh->cubNp + n];
      dfloat c2 = mesh->c2[n + e*mesh->cubNp];

      cubrhsu[n] = -(sigmax-sigmay)*u[n] - sigmay*(sigmay-sigmax)*utilde[n]; // uhat
      cubrhsv[n] = -(sigmay-sigmax)*v[n] - sigmax*(sigmax-sigmay)*vtilde[n]; // vhat
      cubrhsp[n] = -(sigmax+sigmay)*p[n]/c2 - sigmax*sigmay*ptilde[n]; // p

      // update for u~,v~, p~
      cubrhsutilde[n] = u[n]-sigmay*utilde[n]; // du~/dt = -sigmay*u~  + uhat
      cubrhsvtilde[n] = v[n]-sigmax*vtilde[n]; // dv~/dt = -sigmax*v~  + vhat
      cubrhsptilde[n] = c2*p[n];                  // dp~/dt = p
    }

    //Project down and store
    for(iint n=0;n<mesh->Np;++n){
      dfloat rhsp = 0.f;
      dfloat rhsu = 0.f;
      dfloat rhsv = 0.f;
      dfloat rhsptilde = 0.f;
      dfloat rhsutilde = 0.f;
      dfloat rhsvtilde = 0.f;
      for (iint i=0;i<mesh->cubNp;++i){
        rhsp += mesh->cubProject[n*mesh->cubNp + i] * cubrhsp[i];
        rhsu += mesh->cubProject[n*mesh->cubNp + i] * cubrhsu[i];
        rhsv += mesh->cubProject[n*mesh->cubNp + i] * cubrhsv[i];
        rhsptilde += mesh->cubProject[n*mesh->cubNp + i] * cubrhsptilde[i];
        rhsutilde += mesh->cubProject[n*mesh->cubNp + i] * cubrhsutilde[i];
        rhsvtilde += mesh->cubProject[n*mesh->cubNp + i] * cubrhsvtilde[i];
      }

      // store acoustics rhs contributions from PML
      iint rhsBase = 3*mesh->Nfields*(e*mesh->Np + n) + mesh->Nfields*mesh->MRABshiftIndex[lev];
      mesh->rhsq[rhsBase+0] += rhsu;
      mesh->rhsq[rhsBase+1] += rhsv;
      mesh->rhsq[rhsBase+2] += rhsp;

      // store PML rhs
      iint pmlrhsBase = 3*mesh->pmlNfields*(n+pmlId*mesh->Np)  + mesh->pmlNfields*mesh->MRABshiftIndex[lev];
      mesh->pmlrhsq[pmlrhsBase+0] = rhsutilde;
      mesh->pmlrhsq[pmlrhsBase+1] = rhsvtilde;
      mesh->pmlrhsq[pmlrhsBase+2] = rhsptilde;
    }
  }

  free(p); free(u); free(v);
  free(cubrhsp); free(cubrhsu); free(cubrhsv);
  free(ptilde); free(utilde); free(vtilde);
  free(cubrhsptilde); free(cubrhsutilde); free(cubrhsvtilde);
}