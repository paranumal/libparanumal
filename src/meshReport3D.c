
void meshReport3D(const char *mess, mesh3D *mesh){

  printf("%s: (Nfields=%d,Np=%d,Nfaces=%d,Nfp=%d,Nvgeo=%d)\n",
	 mess, mesh->Nfields, mesh->Np, mesh->Nfaces, mesh->Nfp, mesh->Nvgeo);
  
  dfloat maxq = 0, minq = 1e9;
  dfloat maxrhsq = 0, minrhsq = 1e9;

  for(int n=0;n<mesh->Np*mesh->Nelements*mesh->Nfields;++n){
    maxq = mymax(maxq, mesh->q[n]);
    minq = mymin(minq, mesh->q[n]);
    maxrhsq = mymax(maxrhsq, mesh->rhsq[n]);
    minrhsq = mymin(minrhsq, mesh->rhsq[n]);

    printf("%g ", mesh->rhsq[n]);
    if((n%mesh->Nfields) == mesh->Nfields-1)
      printf("\n");
    
  }
  printf("q in %g,%g\n", minq, maxq);
  printf("rhsq in %g,%g\n", minrhsq, maxrhsq);
}
