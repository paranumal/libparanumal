#include "acoustics.h"

dfloat acousticsDopriEstimate(acoustics_t *acoustics){
  
  mesh_t *mesh = acoustics->mesh;
  
  // should insert err = acousticsDopriEstimate2D() here
  //Error estimation 
  //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
  //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
  dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
  acoustics->rkErrorEstimateKernel(Ntotal, 
				   acoustics->ATOL,
				   acoustics->RTOL,
				   acoustics->o_q,
				   acoustics->o_rkq,
				   acoustics->o_rkerr,
				   acoustics->o_errtmp);
  
  acoustics->o_errtmp.copyTo(acoustics->errtmp);
  dfloat localerr = 0;
  dfloat err = 0;
  for(dlong n=0;n<acoustics->Nblock;++n){
    localerr += acoustics->errtmp[n];
  }
  MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

  err = sqrt(err/(mesh->Np*acoustics->totalElements*acoustics->Nfields));
  
  return err;
}
  
