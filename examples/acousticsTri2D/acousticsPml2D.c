#include "mesh2D.h"

// function to compute collocation differentiation
// contributions to nodal DG rhs for acoustics
void acousticsPml2D(mesh2D *mesh){

  // for all elements
  for(iint m=0;m<mesh->pmlNelements;++m){
    iint e = mesh->pmlElementList[m];

    // prefetch geometric factors (constant on triangle)
    dfloat sigmax = mesh->pmlSigmaX[m];
    dfloat sigmay = mesh->pmlSigmaY[m];

    // for all nodes in this element
    for(iint n=0;n<mesh->Np;++n){

      // load state and rhs values at node
      iint   base = mesh->Nfields*(n+e*mesh->Np);
      dfloat u = mesh->q[base+0];
      dfloat v = mesh->q[base+1];
      dfloat p = mesh->q[base+2];

      // load pml state at node
      iint   pmlBase = mesh->pmlNfields*(n+m*mesh->Np);
      dfloat utilde = mesh->pmlq[pmlBase+0];
      dfloat vtilde = mesh->pmlq[pmlBase+1];
      dfloat ptilde = mesh->pmlq[pmlBase+2];

      // update for u,v,p
      dfloat rhsu = -(sigmax-sigmay)*u - sigmay*(sigmay-sigmax)*utilde; // uhat
      dfloat rhsv = -(sigmay-sigmax)*v - sigmax*(sigmax-sigmay)*vtilde; // vhat
      dfloat rhsp = -(sigmax+sigmay)*p - sigmax*sigmay*ptilde; // p

      // update for u~,v~, p~
      dfloat rhsutilde = u-sigmay*utilde; // du~/dt = -sigmay*u~  + uhat
      dfloat rhsvtilde = v-sigmax*vtilde; // dv~/dt = -sigmax*v~  + vhat
      dfloat rhsptilde = p;                  // dp~/dt = p

      // store acoustics rhs contributions from PML
      mesh->rhsq[base+0] += rhsu;
      mesh->rhsq[base+1] += rhsv;
      mesh->rhsq[base+2] += rhsp;

      // store PML rhs
      mesh->pmlrhsq[pmlBase+0] = rhsutilde;
      mesh->pmlrhsq[pmlBase+1] = rhsvtilde;
      mesh->pmlrhsq[pmlBase+2] = rhsptilde;

    }
  }
}
