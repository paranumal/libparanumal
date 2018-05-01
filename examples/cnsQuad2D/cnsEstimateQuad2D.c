#include "cnsQuad2D.h"

dfloat cnsDopriEstimateQuad2D(cns_t *cns){
  
  mesh_t *mesh = cns->mesh;
  
  // should insert err = cnsDopriEstimate2D() here
  //Error estimation 
  //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
  //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
  dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
  cns->rkErrorEstimateKernel(Ntotal, 
			     cns->ATOL,
			     cns->RTOL,
			     cns->o_q,
			     cns->o_rkq,
			     cns->o_rkerr,
			     cns->o_errtmp);
  
  cns->o_errtmp.copyTo(cns->errtmp);
  dfloat localerr = 0;
  dfloat err = 0;
  for(dlong n=0;n<cns->Nblock;++n){
    localerr += cns->errtmp[n];
  }
  MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  // bad normalization
  err = sqrt(err/(cns->totalElements*mesh->Np*mesh->Nfields));
  
  return err;
}
  
