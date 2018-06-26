#include "bns.h"



void bnsRestartWrite(bns_t *bns, dfloat time){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = bns->mesh; 
  // Copy Field To Host
  bns->o_q.copyTo(bns->q);

  // Create Binary File Name
  char fname[BUFSIZ];
  sprintf(fname, "bnsRestartFile_%04d.dat", rank);
  FILE *fp; 
  
  // Overwrite  Binary File
  fp = fopen(fname,"wb");

  // Write q1....qN and all pml variables
  dfloat *elmField  = (dfloat *) calloc(bns->Nfields, sizeof(dfloat)); 
  dfloat *pmlField  = (dfloat *) calloc(bns->Nfields*bns->dim, sizeof(dfloat));

  // First Write the Solution Time
  fwrite(&time, sizeof(dfloat), 1, fp);

  // Write only q works check for MRAB, write history
  for(dlong e = 0; e<mesh->Nelements; e++){
    for(int n=0; n<mesh->Np; n++ ){
      const dlong id = e*bns->Nfields*mesh->Np + n; 
      for(int fld=0; fld<bns->Nfields; fld++){
        elmField[fld] =  bns->q[id + fld*mesh->Np];
      }
      fwrite(elmField, sizeof(dfloat), bns->Nfields, fp);
    }
  }

  if(bns->pmlFlag){
    // Copy Field To Host
    bns->o_pmlqx.copyTo(bns->pmlqx);
    bns->o_pmlqy.copyTo(bns->pmlqy);
    if(bns->dim==3)
      bns->o_pmlqz.copyTo(bns->pmlqz);

    // Write node data
    for(dlong es = 0; es<mesh->pmlNelements; es++)
    {
      dlong e      = mesh->pmlElementIds[es]; 
      dlong pmlId  = mesh->pmlIds[es]; 
      for(int n=0; n<mesh->Np; n++ ){
        const dlong id = pmlId*bns->Nfields*mesh->Np + n;
        for(int fld=0; fld<bns->Nfields; fld++){
          pmlField[fld + 0*bns->Nfields] =  bns->pmlqx[id + fld*mesh->Np];
          pmlField[fld + 1*bns->Nfields] =  bns->pmlqy[id + fld*mesh->Np];
          if(bns->dim==3)
            pmlField[fld + 2*bns->Nfields] =  bns->pmlqz[id + fld*mesh->Np];
        }
      fwrite(pmlField, sizeof(dfloat), bns->Nfields*bns->dim, fp);
      }
    }
  }

  fclose(fp); 
}




void bnsRestartRead(bns_t *bns, dfloat time){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = bns->mesh; 
  // Create Binary File Name
  char fname[BUFSIZ];
  sprintf(fname, "bnsRestartFile_%04d.dat", rank);
  FILE *fp; 
// Keep Writing for PML Variables

  // Overwrite  Binary File
  fp = fopen(fname,"rb");
  
  // Write q1....qN and all pml variables
  dfloat *elmField  = (dfloat *) calloc(bns->Nfields, sizeof(dfloat)); 
  dfloat *pmlField  = (dfloat *) calloc(bns->Nfields*bns->dim, sizeof(dfloat));
  
  dfloat startTime = 0.0; 
  // Update Start Time
  fread(&startTime, sizeof(dfloat),1, fp);

  // Write only q works check for MRAB, write history
  for(dlong e = 0; e<mesh->Nelements; e++){
    for(int n=0; n<mesh->Np; n++ ){
      
      const dlong id = e*bns->Nfields*mesh->Np + n;
      fread(elmField, sizeof(dfloat), bns->Nfields, fp);
      
      for(int fld=0; fld<bns->Nfields; fld++){
        bns->q[id + fld*mesh->Np] = elmField[fld];
      }

    }
  }

  if(bns->pmlFlag){
    
    // Write node data
    for(dlong es = 0; es<mesh->pmlNelements; es++)
    {
      dlong e      = mesh->pmlElementIds[es]; 
      dlong pmlId  = mesh->pmlIds[es]; 
      for(int n=0; n<mesh->Np; n++ ){
        const dlong id = pmlId*bns->Nfields*mesh->Np + n;

        fread(pmlField, sizeof(dfloat), bns->Nfields*bns->dim, fp);
        
        for(int fld=0; fld<bns->Nfields; fld++){          
          bns->pmlqx[id + fld*mesh->Np] = pmlField[fld + 0*bns->Nfields];
          bns->pmlqy[id + fld*mesh->Np] = pmlField[fld + 1*bns->Nfields];
          if(bns->dim==3)
            bns->pmlqz[id + fld*mesh->Np] = pmlField[fld + 2*bns->Nfields];
        }
      }
    } 
  }
 
  

  // Update Fields
  bns->o_q.copyFrom(bns->q);

  if(bns->pmlFlag){
  // Copy Field To Host
    bns->o_pmlqx.copyFrom(bns->pmlqx);
    bns->o_pmlqy.copyFrom(bns->pmlqy);
    if(bns->dim==3)
      bns->o_pmlqz.copyFrom(bns->pmlqz);    
  }

fclose(fp);

// Just Update Time Step Size
bns->startTime = startTime; 
// 











}