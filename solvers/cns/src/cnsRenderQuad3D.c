/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "cns.h"

// interpolate data to plot nodes and save to file (one per process
void cnsRenderQuad3D(cns_t *cns, char *fileBaseName, int fileIndex){

  mesh_t *mesh = cns->mesh;

  int plotNelements = mesh->Nelements*mesh->plotNelements;
  
  dfloat *plotx = (dfloat*) calloc(plotNelements*mesh->Nverts, sizeof(dfloat));
  dfloat *ploty = (dfloat*) calloc(plotNelements*mesh->Nverts, sizeof(dfloat));
  dfloat *plotz = (dfloat*) calloc(plotNelements*mesh->Nverts, sizeof(dfloat));
  dfloat *plotq = (dfloat*) calloc(plotNelements*mesh->Nverts, sizeof(dfloat));
  dfloat *plotVortMag = (dfloat*) calloc(plotNelements*mesh->Nverts, sizeof(dfloat));

  dfloat *tmpPlotx = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  dfloat *tmpPloty = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  dfloat *tmpPlotz = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  dfloat *tmpPlotq = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  dfloat *tmpPlotVortx = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  dfloat *tmpPlotVorty = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  dfloat *tmpPlotVortz = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  dfloat *tmpPlotVortMag = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  
  // compute plot node coordinates on the fly
  int cnt = 0;
  dfloat maxPlotq = -1e9;
  dfloat minPlotq =  1e9;
  for(dlong e=0;e<mesh->Nelements;++e){

    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotxn = 0, plotyn = 0, plotzn=0, plotqn = 0;
      dfloat plotVortxn = 0, plotVortyn = 0, plotVortzn = 0, plotVortMagn = 0;

      for(int m=0;m<mesh->Np;++m){
	dfloat Inm = mesh->plotInterp[n*mesh->Np+m];
	int idm = m+e*mesh->Np;
        plotxn += Inm*mesh->x[idm];
        plotyn += Inm*mesh->y[idm];
        plotzn += Inm*mesh->z[idm];
	plotqn += Inm*cns->q[m+e*mesh->Np*mesh->Nfields];
	plotVortxn += Inm*cns->Vort[m+e*mesh->Np*3 + 0*mesh->Np] ;
	plotVortyn += Inm*cns->Vort[m+e*mesh->Np*3 + 1*mesh->Np] ;
	plotVortzn += Inm*cns->Vort[m+e*mesh->Np*3 + 2*mesh->Np] ;
      }

      maxPlotq = mymax(maxPlotq, plotVortMagn);
      minPlotq = mymin(minPlotq, plotVortMagn);

      tmpPlotx[n] = plotxn;
      tmpPloty[n] = plotyn;
      tmpPlotz[n] = plotzn;
      tmpPlotq[n] = plotqn;
      tmpPlotVortx[n] = plotVortxn;
      tmpPlotVorty[n] = plotVortyn;
      tmpPlotVortz[n] = plotVortzn;

      // mag after interp
#if 0
      tmpPlotVortMag[n] =
	sqrt(plotVortxn*plotVortxn +
	     plotVortyn*plotVortyn +
	     plotVortzn*plotVortzn);
#else
      tmpPlotVortMag[n] =
	(plotxn*plotVortxn +
	 plotyn*plotVortyn +
	 plotzn*plotVortzn);
#endif
    }

    for(int n=0;n<mesh->plotNelements;++n){
      int v1 = mesh->plotEToV[n*mesh->plotNverts+0];
      int v2 = mesh->plotEToV[n*mesh->plotNverts+1];
      int v3 = mesh->plotEToV[n*mesh->plotNverts+2];

      plotx[cnt] = tmpPlotx[v1];
      ploty[cnt] = tmpPloty[v1];
      plotz[cnt] = tmpPlotz[v1];
      plotq[cnt] = tmpPlotq[v1];
      plotVortMag[cnt] = tmpPlotVortMag[v1];
      ++cnt;

      plotx[cnt] = tmpPlotx[v2];
      ploty[cnt] = tmpPloty[v2];
      plotz[cnt] = tmpPlotz[v2];
      plotq[cnt] = tmpPlotq[v2];
      plotVortMag[cnt] = tmpPlotVortMag[v2];
      ++cnt;

      plotx[cnt] = tmpPlotx[v3];
      ploty[cnt] = tmpPloty[v3];
      plotz[cnt] = tmpPlotz[v3];
      plotq[cnt] = tmpPlotq[v3];
      plotVortMag[cnt] = tmpPlotVortMag[v3];
      ++cnt;
      
    }
  }

  int *allNelements = NULL;

  if(mesh->rank==0){
    allNelements = (int*) calloc(mesh->size, sizeof(int));
  }
  
  MPI_Gather(&plotNelements, 1, MPI_INT,
	     allNelements,  1, MPI_INT,
	     0, mesh->comm);
  
  dfloat *allPlotx = NULL;
  dfloat *allPloty = NULL;
  dfloat *allPlotz = NULL;
  dfloat *allPlotq = NULL;
  dfloat *allPlotVortMag = NULL;
  int    *allDispls = NULL;
  
  if(mesh->rank==0){
    
    allDispls = (int*) calloc(mesh->size+1, sizeof(int));
    
    for(int r=0;r<mesh->size;++r){
      allNelements[r] *= 3;
    }
    for(int r=1;r<=mesh->size;++r){
      allDispls[r] = allDispls[r-1]+allNelements[r-1];
    }

    allPlotx  = (dfloat*) calloc(allDispls[mesh->size], sizeof(dfloat));
    allPloty  = (dfloat*) calloc(allDispls[mesh->size], sizeof(dfloat));
    allPlotz  = (dfloat*) calloc(allDispls[mesh->size], sizeof(dfloat));
    allPlotq  = (dfloat*) calloc(allDispls[mesh->size], sizeof(dfloat));
    allPlotVortMag  = (dfloat*) calloc(allDispls[mesh->size], sizeof(dfloat));
  }
  
  MPI_Gatherv(plotx, plotNelements*3, MPI_DFLOAT,
	      allPlotx, allNelements, allDispls, MPI_DFLOAT, 0, mesh->comm);

  MPI_Gatherv(ploty, plotNelements*3, MPI_DFLOAT,
	      allPloty, allNelements, allDispls, MPI_DFLOAT, 0, mesh->comm);

  MPI_Gatherv(plotz, plotNelements*3, MPI_DFLOAT,
	      allPlotz, allNelements, allDispls, MPI_DFLOAT, 0, mesh->comm);

  MPI_Gatherv(plotq, plotNelements*3, MPI_DFLOAT,
	      allPlotq, allNelements, allDispls, MPI_DFLOAT, 0, mesh->comm);

  MPI_Gatherv(plotVortMag, plotNelements*3, MPI_DFLOAT,
	      allPlotVortMag, allNelements, allDispls, MPI_DFLOAT, 0, mesh->comm);

  if(mesh->rank==0){
    simpleRayTracer(allDispls[mesh->size]/3,
		    allPlotx,
		    allPloty,
		    allPlotz,
		    allPlotVortMag,
		    fileBaseName,
		    fileIndex);
    
    printf("plotq in [%g,%g]\n", minPlotq, maxPlotq);
    free(allPlotx);
    free(allPloty);
    free(allPlotz);
    free(allPlotq);
    free(allPlotVortMag);
    free(allDispls);
  }
  
  free(tmpPlotx);
  free(tmpPloty);
  free(tmpPlotz);
  free(tmpPlotVortx);
  free(tmpPlotVorty);
  free(tmpPlotVortz);
  free(tmpPlotq);
  free(tmpPlotVortMag);

  free(plotx);
  free(ploty);
  free(plotz);
  free(plotq);
  free(plotVortMag);
}
