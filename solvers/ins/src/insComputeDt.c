#include "ins.h"

void insComputeDt(ins_t *ins, dfloat time){

  mesh_t *mesh = ins->mesh; 
  // copy data to host
  ins->o_U.copyTo(ins->U);

  dfloat hminL = 0.0, umaxL = 0.0, dt = 1e9;
  for(dlong e=0;e<mesh->Nelements;++e){
    hminL = 1e9; 
    umaxL = 0.0; 
    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      dfloat hest = 2./(sJ*invJ);
      hminL = mymin(hminL, hest);
    }

    for(int n=0;n<mesh->Np;++n){
      const dlong id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat uxn = ins->U[id+0*ins->fieldOffset];
      dfloat uyn = ins->U[id+1*ins->fieldOffset];
      dfloat uzn = 0.0;
      if (ins->dim==3) 
        uzn = ins->U[id+2*ins->fieldOffset];

      //Squared maximum velocity
      dfloat numax;
      if (ins->dim==2)
        numax = uxn*uxn + uyn*uyn;
      else 
        numax = uxn*uxn + uyn*uyn + uzn*uzn;

      umaxL = mymax(umaxL, numax);
    }
    umaxL = sqrt(umaxL);
    
    //Guard for around zero velocity
    umaxL = (umaxL<1.E-12) ? 1.E-3 : umaxL;

    dfloat dtL = ins->cfl*hminL/((mesh->N+1)*(mesh->N+1)*umaxL);
    dt = mymin(dt, dtL); 
  }


// Save the time step size
// ins->dto = ins->dt; 
// MPI_Allreduce to get global minimum dt
MPI_Allreduce(&dt, &(ins->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

// Update dt dependent variables 
ins->idt    = 1.0/ins->dt;
ins->lambda = ins->g0 / (ins->dt * ins->nu);

}