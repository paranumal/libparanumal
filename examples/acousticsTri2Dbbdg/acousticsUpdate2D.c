#include "mesh2D.h"

void acousticsUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb){
  
  // Low storage Runge Kutta time step update
  for(iint n=0;n<mesh->Nelements*mesh->Np*mesh->Nfields;++n){

    mesh->resq[n] = rka*mesh->resq[n] + mesh->dt*mesh->rhsq[n];

    mesh->q[n] += rkb*mesh->resq[n];
  }
}

void acousticsUpdate2D_wadg(mesh2D *mesh,
                            dfloat rka,
                            dfloat rkb){

  // Low storage Runge Kutta time step update
  for(iint e=0;e<mesh->Nelements;++e){  // for all elements
    dfloat p[mesh->cubNp];

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp;++n){
      p[n] = 0.f;
      iint id = mesh->Nfields*(e*mesh->Np);
      for (iint i=0;i<mesh->Np;++i){
        p[n] += mesh->cubInterp[n*mesh->Np + i] * mesh->rhsq[id + mesh->Nfields*i + 2];
      }

      // Multiply result by wavespeed c2 at cubature node
      p[n] *= mesh->c2[n + e*mesh->cubNp];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np;++n){
      // Extract velocity rhs
      iint id = mesh->Nfields*(e*mesh->Np + n);
      dfloat rhsqn[mesh->Nfields];
      rhsqn[0] = mesh->rhsq[id + 0];
      rhsqn[1] = mesh->rhsq[id + 1];  
      
      // Project scaled rhs down
      dfloat rhsp = 0.f;
      for (iint i=0;i<mesh->cubNp;++i){
        rhsp += mesh->cubProject[n*mesh->cubNp + i] * p[i];
      }
      //rhsqn[2] = mesh->rhsq[id + 2];
      rhsqn[2] = rhsp;

      // Increment solutions
      for (iint fld = 0; fld < mesh->Nfields; ++fld){ 
        mesh->resq[id + fld] = rka*mesh->resq[id + fld] + mesh->dt*rhsqn[fld];
        mesh->q[id + fld] += rkb*mesh->resq[id + fld];
      }      
    }
  }
}