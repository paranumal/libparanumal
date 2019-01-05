#if 0
  dfloat *q1 = (dfloat*) calloc(mesh->Np*mesh->Nelements,sizeof(dfloat));
  dfloat *q2 = (dfloat*) calloc(mesh->Np*mesh->Nelements,sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Nelements;++n){
    q1[n] = pow(mesh->x[n],9); 
    q2[n] = q1[n];
  }
  occa::memory o_q1   = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), q1);
  occa::memory o_q2   = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), q2);
  occa::memory o_Aq1 = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), q1);
  occa::memory o_Aq2 = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), q2);

  elliptic->partialAxKernel(mesh->NlocalGatherElements,
			    mesh->o_localGatherElementList,
			    mesh->o_ggeo,
			    mesh->o_Dmatrices,
			    mesh->o_Smatrices,
			    mesh->o_MM,
			    lambda,
			    o_q1,
			    o_Aq1);
  
  elliptic->partialCubatureAxKernel(mesh->NlocalGatherElements,
				    mesh->o_localGatherElementList,
				    mesh->o_cubggeo,
				    mesh->o_cubD,
				    mesh->o_cubInterpT,
				    lambda, o_q2, o_Aq2);


  o_Aq1.copyTo(q1);
  o_Aq2.copyTo(q2);

  for(int n=0;n<mesh->Np*mesh->Nelements;++n){
    printf("Aq1[%d] = %lf, Aq2[%d] = %lf, diff[%d] = %lg\n",
	   n, q1[n], n, q2[n], n, fabs(q1[n]-q2[n]));
  }

  
  exit(-1);
#endif
  
