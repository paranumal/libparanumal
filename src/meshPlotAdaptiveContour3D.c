#include "mesh3D.h"

void PlotAdaptiveContour3D(mesh_t *mesh, char *fname, dfloat *u, int Nlevels, dfloat *levels, dfloat tol){

  // function PlotAdaptiveContour3D(u, levels, tol)
  // Purpose: adaptively refine the mesh to approximately locate isocontours


  // build interpolation matrix (coarse->fine)
  // assume these are loaded from node file
  // mesh->contourEToV = [1 5 7 8; 5 2 6 9; 7 6 3 10; 8 9 10 4; 8 5 7 9; 7 5 6 9; 8 9 7 10; 9 6 7 10];
  // mesh->contourVX   = [-1  1 -1 -1  0  0 -1 -1  0 -1];
  // mesh->contourVY   = [-1 -1  1 -1 -1  0  0 -1 -1  0];
  // mesh->contourVZ   = [-1 -1 -1  1 -1 -1 -1  0  0  0];
  // mesh->contourInterpN
  // v1 = EToVi(:,1); v2 = EToVi(:,2); v3 = EToVi(:,3); v4 = EToVi(:,4);
  // ri = 0.5*(-(r+s+t+1)*VXi(v1) + (1+r)*VXi(v2) + (1+s)*VXi(v3) + (1+t)*VXi(v4) );
  // si = 0.5*(-(r+s+t+1)*VYi(v1) + (1+r)*VYi(v2) + (1+s)*VYi(v3) + (1+t)*VYi(v4) );
  // ti = 0.5*(-(r+s+t+1)*VZi(v1) + (1+r)*VZi(v2) + (1+s)*VZi(v3) + (1+t)*VZi(v4) );
  //interp = Vandermonde3D(N, ri(:), si(:), ti(:))*invV;

  // mesh->contourInterp1
  // ri = [-1;1;-1;-1]; si = [-1;-1;1;-1]; ti = [-1;-1;-1;1]; refNp = length(ri);
  // interp1 = Vandermonde3D(N, ri(:), si(:), ti(:))*invV;
  // mesh->contourF
  //sk = 1;
  //F = spalloc(Np,Np,1);
  //for i=0:N % old ordering
  //  for j=0:N - i
  //    for k=0:N - i - j
  //      if(i+j+k<=1), F(sk,sk) = 1.; end;
  //      sk = sk+1;
  //    end
  //  end
  //end

  // contourFilter:     ufilt = V*F*invV

  int totalNelements = 0;
  dfloat *plotx = (dfloat*) calloc(4, sizeof(dfloat));
  dfloat *ploty = (dfloat*) calloc(4, sizeof(dfloat));
  dfloat *plotz = (dfloat*) calloc(4, sizeof(dfloat));
  dfloat *plotu = (dfloat*) calloc(4, sizeof(dfloat));
  int plotNp = 4;

  for(int lev=1;lev<=Nlevels;++lev){
    
    int Nelements = mesh->Nelements;
    int Np = mesh->Np;
    
    dfloat *refu = (dfloat*) calloc(Nelements*Np, sizeof(dfloat));
    dfloat *refx = (dfloat*) calloc(Nelements*Np, sizeof(dfloat));
    dfloat *refy = (dfloat*) calloc(Nelements*Np, sizeof(dfloat));
    dfloat *refz = (dfloat*) calloc(Nelements*Np, sizeof(dfloat));
    
    for(int n=0;n<Np*Nelements;++n){
      refu[n] = u[n];
      refx[n] = mesh->x[n];
      refy[n] = mesh->y[n];
      refz[n] = mesh->z[n];
    }
    
    dfloat *newu, *newx, *newy, *newz, *newJ;
    
    dfloat err = 1;
    while(err > tol){ // should add max refinement check here
      
      dfloat level = levels[lev-1];
      
      int *refineList = (int*) calloc(Nelements,sizeof(int));
      int Nrefine = 0;
      for(int e=0;e<Nelements;++e){
	dfloat umin = refu[e*Np+0];
	dfloat umax = refu[e*Np+0];
	
	for(int n=1;n<Np;++n){
	  umin = mymin(umin, refu[e*Np+n]);
	  umax = mymax(umax, refu[e*Np+n]);
	}
	
	if(umin<=level && umax>=level){
	  refineList[Nrefine] = e;
	  ++Nrefine;
	}
      }
      
      int newNelements = 8*Nrefine;

      newu = (dfloat*) calloc(Np*newNelements, sizeof(dfloat));
      newx = (dfloat*) calloc(Np*newNelements, sizeof(dfloat));
      newy = (dfloat*) calloc(Np*newNelements, sizeof(dfloat));
      newz = (dfloat*) calloc(Np*newNelements, sizeof(dfloat));
      newJ = (dfloat*) calloc(Np*newNelements, sizeof(dfloat));

      for(int n=0;n<Nrefine;++n){
	int e = refineList[n];
	for(int m=0;m<8*Np;++m){
	  for(int i=0;i<Np;++i){
	    // note layout
	    newu[8*Np*n+m] += mesh->contourInterp[m*Np + i]*refu[e*Np+i];
	    newx[8*Np*n+m] += mesh->contourInterp[m*Np + i]*refx[e*Np+i];
	    newy[8*Np*n+m] += mesh->contourInterp[m*Np + i]*refy[e*Np+i];
	    newz[8*Np*n+m] += mesh->contourInterp[m*Np + i]*refz[e*Np+i];
	  }
	}
      }
      
      free(refu);
      free(refx);
      free(refy);
      free(refz);

      Nelements = newNelements;
      refu = newu;
      refx = newx;
      refy = newy;
      refz = newz;

      err = 0;
      for(int e=0;e<Nelements;++e){
	for(int n=0;n<Np;++n){
	  dfloat errn = -refu[e*Np+n];
	  for(int m=0;m<Np;++m)
	    errn += mesh->contourFilter[n*Np+m]*refu[e*Np+m];
	  err = mymax(err, fabs(errn));
	}
      }
    }
    
    // append to lists
    plotx = (dfloat*) realloc(plotx, 4*(totalNelements+Nelements)*sizeof(dfloat));
    ploty = (dfloat*) realloc(ploty, 4*(totalNelements+Nelements)*sizeof(dfloat));
    plotz = (dfloat*) realloc(plotz, 4*(totalNelements+Nelements)*sizeof(dfloat));
    plotu = (dfloat*) realloc(plotu, 4*(totalNelements+Nelements)*sizeof(dfloat));
    
    for(int e=0;e<Nelements;++e){
      for(int n=0;n<plotNp;++n){
	
	dfloat px = 0, py = 0, pz = 0, pu = 0;
	
	for(int m=0;m<Np;++m){
	  px += mesh->contourInterp1[n*Np+m]*refx[e*Np+m];
	  py += mesh->contourInterp1[n*Np+m]*refy[e*Np+m];
	  pz += mesh->contourInterp1[n*Np+m]*refz[e*Np+m];
	  pu += mesh->contourInterp1[n*Np+m]*refu[e*Np+m];
	}
	
	plotx[(e+totalNelements)*plotNp+n] = px;
	ploty[(e+totalNelements)*plotNp+n] = py;
	plotz[(e+totalNelements)*plotNp+n] = pz;
	plotu[(e+totalNelements)*plotNp+n] = pu;
	
      }
    }
    
    totalNelements += Nelements;
    
    free(refx);
    free(refy);
    free(refz);
    free(refu);
  }
  
  int plotNelements = totalNelements;

  FILE *fp = fopen(fname, "w");
  
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
	  plotNelements*plotNp,
	  plotNelements);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(iint n=0;n<plotNelements*plotNp;++n){
    fprintf(fp, "       ");
    fprintf(fp, "%g %g %g\n", plotx[n],ploty[n],plotz[n]);
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");
  
  // write out pressure
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"pressure\" Format=\"ascii\">\n");
  
  for(iint e=0;e<plotNelements;++e){
    for(iint n=0;n<plotNp;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotu[e*plotNp+n]);
    }
  }
  
  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(iint e=0;e<plotNelements;++e){
    fprintf(fp, "       ");
    for(int m=0;m<plotNverts;++m){
      fprintf(fp, "%d ", e*plotNp + m);
    }
    fprintf(fp, "\n");
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  iint cnt = 0;
  for(iint e=0;e<plotNelements;++e){
    cnt += plotNverts;
    fprintf(fp, "       ");
    fprintf(fp, "%d\n", cnt);
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(iint e=0;e<plotNelements;++e){
    fprintf(fp, "10\n"); // TET code ?
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
  
  fclose(fp);
  
}
