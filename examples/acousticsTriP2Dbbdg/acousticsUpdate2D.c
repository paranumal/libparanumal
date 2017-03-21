#include "acoustics2D.h"

void acousticsUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb){
  
	// Low storage Runge Kutta time step update
	for(iint n=0;n<mesh->Nelements*mesh->NpMax*mesh->Nfields;++n){

		mesh->resq[n] = rka*mesh->resq[n] + mesh->dt*mesh->rhsq[n];

		mesh->q[n] += rkb*mesh->resq[n];
	}
}

void acousticsUpdate2D_wadg(mesh2D *mesh,
                            dfloat rka,
                            dfloat rkb){

  // Low storage Runge Kutta time step update
  for(iint e=0;e<mesh->Nelements;++e){  // for all elements
    iint N = mesh->N[e];

    dfloat p[mesh->cubNpMax];

    // Interpolate rhs to cubature nodes
    for(iint n=0;n<mesh->cubNp[N];++n){
      p[n] = 0.f;
      iint id = mesh->Nfields*(e*mesh->NpMax);
      for (iint i=0;i<mesh->Np[N];++i){
        p[n] += mesh->cubInterp[N][n*mesh->Np[N] + i] * mesh->rhsq[id + mesh->Nfields*i + 2];
      }

      // Multiply result by wavespeed c2 at cubature node
      p[n] *= mesh->c2[n + e*mesh->cubNpMax];
    }

    // Increment solution, project result back down
    for(iint n=0;n<mesh->Np[N];++n){
      // Extract velocity rhs
      iint id = mesh->Nfields*(e*mesh->NpMax + n);
      dfloat rhsqn[mesh->Nfields];
      rhsqn[0] = mesh->rhsq[id + 0];
      rhsqn[1] = mesh->rhsq[id + 1];  
      
      // Project scaled rhs down
      dfloat rhsp = 0.f;
      for (iint i=0;i<mesh->cubNp[N];++i){
        rhsp += mesh->cubProject[N][n*mesh->cubNp[N] + i] * p[i];
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