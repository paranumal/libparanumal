#include <stdlib.h>
#include "mesh.h"
#include <math.h>

//original mapping (r,s) -> (x,y)
void meshPhysicalNodesQuad3D(mesh_t *mesh){
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];

    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];

    
    for(iint n=0;n<mesh->Np;++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* physical coordinate of interpolation node */
      dfloat xlin = 
	+0.25*(1-rn)*(1-sn)*xe1
	+0.25*(1+rn)*(1-sn)*xe2
	+0.25*(1+rn)*(1+sn)*xe3
	+0.25*(1-rn)*(1+sn)*xe4;

      dfloat ylin =
	+0.25*(1-rn)*(1-sn)*ye1
	+0.25*(1+rn)*(1-sn)*ye2
	+0.25*(1+rn)*(1+sn)*ye3
	+0.25*(1-rn)*(1+sn)*ye4;

      dfloat zlin =
	+0.25*(1-rn)*(1-sn)*ze1
	+0.25*(1+rn)*(1-sn)*ze2
	+0.25*(1+rn)*(1+sn)*ze3
	+0.25*(1-rn)*(1+sn)*ze4;

      //      printf("xlin=%g, ylin=%g, zlin=%g\n", xlin, ylin, zlin);
      
      // project to sphere
      dfloat rlin = sqrt(xlin*xlin+ylin*ylin+zlin*zlin);
      mesh->x[cnt] = xlin/rlin; 
      mesh->y[cnt] = ylin/rlin; 
      mesh->z[cnt] = zlin/rlin; 
	
      ++cnt;
    }
  }
}

//original mapping (r,s) -> (x,y)
void meshSphericalNodesQuad3D(mesh_t *mesh){

  //constants used in conversions
  const dfloat R = 1;
  const dfloat a = 1./sqrt(3.);
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];

    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];

    int faceNumber = mesh->cubeFaceNumber[e];
    
    for(iint n=0;n<mesh->Np;++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* physical coordinate of node on cube */
      dfloat xlin = 
	+0.25*(1-rn)*(1-sn)*xe1
	+0.25*(1+rn)*(1-sn)*xe2
	+0.25*(1+rn)*(1+sn)*xe3
	+0.25*(1-rn)*(1+sn)*xe4;

      dfloat ylin =
	+0.25*(1-rn)*(1-sn)*ye1
	+0.25*(1+rn)*(1-sn)*ye2
	+0.25*(1+rn)*(1+sn)*ye3
	+0.25*(1-rn)*(1+sn)*ye4;

      dfloat zlin =
	+0.25*(1-rn)*(1-sn)*ze1
	+0.25*(1+rn)*(1-sn)*ze2
	+0.25*(1+rn)*(1+sn)*ze3
	+0.25*(1-rn)*(1+sn)*ze4;

      dfloat xsph, ysph, zsph, cubRad;

      cubRad = sqrt(xlin*xlin + ylin*ylin + zlin*zlin);
      xsph = R/cubRad * xlin;
      ysph = R/cubRad * ylin;
      zsph = R/cubRad * zlin;

      if (fabs(xsph*xsph + ysph*ysph + zsph*zsph - 1) > 1e-12) printf("error\n");
	  
      // project to sphere
      mesh->x[cnt] = zsph; 
      mesh->y[cnt] = ysph; 
      mesh->z[cnt] = -xsph;
      
      ++cnt;
    }
  }
}

//original mapping (r,s) -> (x,y)
void meshEquiSphericalNodesQuad3D(mesh_t *mesh){

  //constants used in conversions
  const dfloat R = 1;
  const dfloat a = 1./sqrt(3.);
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

  mesh->rphysical = (dfloat*) calloc(mesh->NgridElements*mesh->Np,sizeof(dfloat));
  mesh->sphysical = (dfloat*) calloc(mesh->NgridElements*mesh->Np,sizeof(dfloat));

  mesh->rlocal = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->slocal = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];

    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];

    int faceNumber = mesh->cubeFaceNumber[e];
    
    for(iint n=0;n<mesh->Np;++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* physical coordinate of node on cube */
      dfloat xlin = 
	+0.25*(1-rn)*(1-sn)*xe1
	+0.25*(1+rn)*(1-sn)*xe2
	+0.25*(1+rn)*(1+sn)*xe3
	+0.25*(1-rn)*(1+sn)*xe4;

      dfloat ylin =
	+0.25*(1-rn)*(1-sn)*ye1
	+0.25*(1+rn)*(1-sn)*ye2
	+0.25*(1+rn)*(1+sn)*ye3
	+0.25*(1-rn)*(1+sn)*ye4;

      dfloat zlin =
	+0.25*(1-rn)*(1-sn)*ze1
	+0.25*(1+rn)*(1-sn)*ze2
	+0.25*(1+rn)*(1+sn)*ze3
	+0.25*(1-rn)*(1+sn)*ze4;

      dfloat xsph, ysph, zsph, norm, rlin, slin;

      dfloat faceScale = sqrt(3)*M_PI/4;
      
      switch(faceNumber) {
      case 1: //positive constant x
	xsph = cos(faceScale*ylin)*cos(faceScale*zlin);
	ysph = sin(faceScale*ylin)*cos(faceScale*zlin);
	zsph = sin(faceScale*zlin)*cos(faceScale*ylin);

	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);

	xsph /= norm;
	ysph /= norm;
	zsph /= norm;
	
	rlin = faceScale*ylin;
        slin = faceScale*zlin;
	
	break;
      case 2: //positive constant y
	xsph = sin(faceScale*xlin)*cos(faceScale*zlin);
	ysph = cos(faceScale*xlin)*cos(faceScale*zlin);
	zsph = sin(faceScale*zlin)*cos(faceScale*xlin);

	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);

	xsph /= norm;
	ysph /= norm;
	zsph /= norm;

	rlin = faceScale*xlin;
	slin = faceScale*zlin;
	
	break;
      case 3: //negative constant x
	xsph = -1*cos(faceScale*ylin)*cos(faceScale*zlin);
	ysph = sin(faceScale*ylin)*cos(faceScale*zlin);
	zsph = sin(faceScale*zlin)*cos(faceScale*ylin);

	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);

	xsph /= norm;
	ysph /= norm;
	zsph /= norm;

	rlin = faceScale*ylin;
	slin = faceScale*zlin;
	
	break;
      case 4: //negative constant y
	xsph = sin(faceScale*xlin)*cos(faceScale*zlin);
	ysph = -1*cos(faceScale*xlin)*cos(faceScale*zlin);
	zsph = sin(faceScale*zlin)*cos(faceScale*xlin);

	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);

	xsph /= norm;
	ysph /= norm;
	zsph /= norm;

	rlin = faceScale*xlin;
	slin = faceScale*zlin;
	
	break;
      case 5: //positive constant z
	xsph = sin(faceScale*xlin)*cos(faceScale*ylin);
	ysph = sin(faceScale*ylin)*cos(faceScale*xlin);
	zsph = cos(faceScale*xlin)*cos(faceScale*ylin);

	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);

	xsph /= norm;
	ysph /= norm;
	zsph /= norm;

	rlin = faceScale*xlin;
	slin = faceScale*ylin;
	
	break;
      case 6: //negative constant z
	xsph =  sin(faceScale*xlin)*cos(faceScale*ylin);
	ysph =  sin(faceScale*ylin)*cos(faceScale*xlin);
	zsph = -1*cos(faceScale*xlin)*cos(faceScale*ylin);

	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);

	xsph /= norm;
	ysph /= norm;
	zsph /= norm;

        rlin = faceScale*xlin;
	slin = faceScale*ylin;
	
	break;
      }

      if (fabs(xsph*xsph + ysph*ysph + zsph*zsph - 1) > 1e-12) printf("error %d\n",mesh->cubeFaceNumber[e]);
	  
      // project to sphere
      mesh->x[cnt] = xsph; 
      mesh->y[cnt] = ysph; 
      mesh->z[cnt] = zsph;

      mesh->rlocal[cnt] = (mesh->r[n]+1)/2;
      mesh->slocal[cnt] = (mesh->s[n]+1)/2;

      mesh->rphysical[cnt] = rlin;
      mesh->sphysical[cnt] = slin;
      
      ++cnt;
    }
    /*if (e == 0) {
      for (iint n = 0; n < mesh->Np; ++n) {
	printf("%d %lf %lf\n",n,mesh->rlocal[e*mesh->Np + n],mesh->slocal[e*mesh->Np + n]);
      }
      }*/

    /*for (iint p = 0; p < mesh->Nq; ++p) {
      for (iint q = 0; q < mesh->Nq; ++q) {
	for (iint t = 0; t < mesh->Nq; ++t) {
	  if ((q != t) && (fabs(mesh->slocal[e*mesh->Np + q*mesh->Nq + p] - mesh->slocal[e*mesh->Np + t*mesh->Nq + p]) < 1e-10)) {
	    printf("slocal should not be equal\n");
	  }
	  if ((q != t) && (fabs(mesh->rlocal[e*mesh->Np + p*mesh->Nq + q] - mesh->rlocal[e*mesh->Np + p*mesh->Nq + t]) < 1e-10)) {
	    printf("rlocal should not be equal\n");
	  }
	}
      }
      }*/
  }
}

void meshEquiSphericalExtensionQuad3D(mesh_t *mesh) {

  mesh->eInterp = (iint *) calloc((mesh->NgridElements - mesh->Nelements)*mesh->Np,sizeof(iint));
  mesh->overlapDirection = (char *) calloc(mesh->NgridElements - mesh->Nelements,sizeof(char));
  mesh->perp_index = (iint *) calloc((mesh->NgridElements - mesh->Nelements)*mesh->Np,sizeof(iint));
  mesh->par_loc = (dfloat *) calloc((mesh->NgridElements - mesh->Nelements)*mesh->Np,sizeof(dfloat));
  
  for (iint i = 0; i < mesh->NgridElements - mesh->Nelements; ++i) {
    iint eOverlap = mesh->overlap[i];
    
    iint eAdj;
    for (iint f = 0; f < mesh->Nfaces; ++f) {
      if (mesh->cubeFaceNumber[eOverlap] != mesh->cubeFaceNumber[mesh->EToE[eOverlap*mesh->Nfaces + f]]) {
	eAdj = mesh->EToE[eOverlap*mesh->Nfaces + f];
	break;
      }
    }
    
    iint e = mesh->Nelements + i;
    
    iint face1 = mesh->cubeFaceNumber[eOverlap] - 1;
    iint face2 = mesh->cubeFaceNumber[eAdj] - 1;
    
    for (iint n = 0; n < mesh->Np; ++n) {
      
      dfloat offset = M_PI/(2*mesh->edgeLength);
      
      iint faceHash = 6*face1 + face2;
      
      dfloat rold, sold, rnew, snew, par_loc, perp_index, rmin, smin, rmax, smax;
      iint rsdir, eInterp;
      iint rdir = 0;
      iint sdir = 1;
      dfloat tol = 1e-10;
      
      switch (faceHash)
	{
	case 13:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] - offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = rold + M_PI/2;
	  snew = -1*atan(tan(sold)/tan(rold));

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  perp_index = (mesh->Np - n - 1)%mesh->Nq;
	  par_loc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;
	  
	  break;
	case 20:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] - offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = -1*(rold + M_PI/2);
	  snew = -1*atan(tan(sold)/tan(rold));

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  perp_index = n%mesh->Nq;
	  par_loc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;
	  
	  break;
	case 6:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] + offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = M_PI/2 - rold;
	  snew = atan(tan(sold)/tan(rold));

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  perp_index = n%mesh->Nq;
	  par_loc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;

	  break;
	case 3:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] + offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = rold - M_PI/2;
	  snew = atan(tan(sold)/tan(rold));

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  perp_index = (mesh->Np - n - 1)%mesh->Nq;
	  par_loc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;

	  break;
	case 1:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] + offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = rold - M_PI/2;
	  snew = atan(tan(sold)/tan(rold));

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  perp_index = n%mesh->Nq;
	  par_loc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;

	  break;
	case 8:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] + offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = rold - M_PI/2;
	  snew = atan(tan(sold)/tan(rold));

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  perp_index = n%mesh->Nq;
	  par_loc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;

	  break;
	case 15:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] - offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = -1*(rold + M_PI/2);
	  snew = -1*atan(tan(sold)/tan(rold));

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];

	  else eInterp = eOverlap;
	  
	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  perp_index = (mesh->Np - n - 1)%mesh->Nq;
	  par_loc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;

	  break;
	case 18:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] - offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = rold + M_PI/2;
	  snew = -1*atan(tan(sold)/tan(rold));

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];

	  else eInterp = eOverlap;
	  
	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  perp_index = n%mesh->Nq;
	  par_loc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;

	  break;
	case 4:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] + offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = atan(tan(sold)/tan(rold));
	  snew = M_PI/2 - rold;

	  rmin = mesh->rphysical[eOverlap*mesh->Np];
	  rmax = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  rmin = mesh->rphysical[eInterp*mesh->Np];
	  rmax = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = (mesh->Np - n -1)/mesh->Nq;
	  
	  rsdir = rdir;

	  break;
	case 30:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] - offset;

	  snew = -1*atan(tan(rold)/tan(sold));
	  rnew = sold + M_PI/2;

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (snew - smin)/(smax - smin);
	  perp_index = (mesh->Np - n - 1)%mesh->Nq;
	  
	  rsdir = sdir;

	  break;
	case 5:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] + offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = atan(tan(sold)/tan(rold));
	  snew = rold - M_PI/2;

	  rmin = mesh->rphysical[eOverlap*mesh->Np];
	  rmax = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  rmin = mesh->rphysical[eInterp*mesh->Np];
	  rmax = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = n/mesh->Nq;
	  
	  rsdir = rdir;

	  break;
	case 24:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] + offset;

	  snew = atan(tan(rold)/tan(sold));
	  rnew = M_PI/2 - sold;

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (snew - smin)/(smax - smin);
	  perp_index = (mesh->Np - n - 1)%mesh->Nq;
	  
	  rsdir = sdir;

	  break;
	case 16:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] - offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = -1*atan(tan(sold)/tan(rold));
	  snew = M_PI/2. + rold;

	  rmax = mesh->rphysical[eOverlap*mesh->Np];
	  rmin = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else eInterp = eOverlap;

	  rmax = mesh->rphysical[eInterp*mesh->Np];
	  rmin = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];

	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = (mesh->Np - n - 1)/mesh->Nq;
	  
	  rsdir = rdir;

	  break;
	case 26:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] + offset;

	  snew = atan(tan(rold)/tan(sold));
	  rnew = sold - M_PI/2;

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (snew - smin)/(smax - smin);
	  perp_index = n%mesh->Nq;
	  
	  rsdir = sdir;

	  break;
	case 17:
	  rold = mesh->rphysical[eAdj*mesh->Np + n] - offset;
	  sold = mesh->sphysical[eAdj*mesh->Np + n];

	  rnew = -1*atan(tan(sold)/tan(rold));
	  snew = -1*(rold + M_PI/2.);

	  rmax = mesh->rphysical[eOverlap*mesh->Np];
	  rmin = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else eInterp = eOverlap;

	  rmax = mesh->rphysical[eInterp*mesh->Np];
	  rmin = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = n/mesh->Nq;
	  
	  rsdir = rdir;

	  break;
	case 32:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] - offset;

	  snew = -1*atan(tan(rold)/tan(sold));
	  rnew = -1*(M_PI/2. + sold);

	  smin = mesh->sphysical[eOverlap*mesh->Np];
	  smax = mesh->sphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  smin = mesh->sphysical[eInterp*mesh->Np];
	  smax = mesh->sphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (snew - smin)/(smax - smin);
	  perp_index = n%mesh->Nq;
	  
	  rsdir = sdir;

	  break;
	case 25:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] + offset;

	  snew = M_PI/2. - sold;
	  rnew = atan(tan(rold)/tan(sold));

	  rmax = mesh->rphysical[eOverlap*mesh->Np];
	  rmin = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else eInterp = eOverlap;

	  rmax = mesh->rphysical[eInterp*mesh->Np];
	  rmin = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = (mesh->Np - n - 1)/mesh->Nq;
	  
	  rsdir = rdir;

	  break;
	case 10:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] + offset;

	  snew = M_PI/2. - sold;
	  rnew = atan(tan(rold)/tan(sold));

	  rmax = mesh->rphysical[eOverlap*mesh->Np];
	  rmin = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else eInterp = eOverlap;
	  
	  rmax = mesh->rphysical[eInterp*mesh->Np];
	  rmin = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = (mesh->Np - n - 1)/mesh->Nq;
	  
	  rsdir = rdir;

	  break;
	case 22:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] - offset;

	  snew = sold + M_PI/2.;
	  rnew = -1*atan(tan(rold)/tan(sold));

	  rmin = mesh->rphysical[eOverlap*mesh->Np];
	  rmax = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  rmin = mesh->rphysical[eInterp*mesh->Np];
	  rmax = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = (mesh->Np - n - 1)/mesh->Nq;
	  
	  rsdir = rdir;

	  break; 
	case 27:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] + offset;

	  snew = sold - M_PI/2.;
	  rnew = atan(tan(rold)/tan(sold));

	  rmax = mesh->rphysical[eOverlap*mesh->Np];
	  rmin = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else eInterp = eOverlap;

	  rmax = mesh->rphysical[eInterp*mesh->Np];
	  rmin = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = n/mesh->Nq;
	  
	  rsdir = rdir;

	  break;
	case 11:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] + offset;

	  snew = sold - M_PI/2.;
	  rnew = atan(tan(rold)/tan(sold));

	  rmax = mesh->rphysical[eOverlap*mesh->Np];
	  rmin = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else eInterp = eOverlap;

	  rmax = mesh->rphysical[eInterp*mesh->Np];
	  rmin = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = n/mesh->Nq;
	  
	  rsdir = rdir;

	  break;
	case 31:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] - offset;

	  snew = M_PI/2. + sold;
	  rnew = -1*atan(tan(rold)/tan(sold));

	  rmin = mesh->rphysical[eOverlap*mesh->Np];
	  rmax = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;
	  
	  rmin = mesh->rphysical[eInterp*mesh->Np];
	  rmax = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = (mesh->Np - n - 1)/mesh->Nq;
	  
	  rsdir = rdir;

	  break;
	case 23:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] - offset;

	  rnew = -1*(sold + M_PI/2.);
	  snew = -1*atan(tan(rold)/tan(sold));

	  smin = mesh->rphysical[eOverlap*mesh->Np];
	  smax = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (smin - snew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (snew - smax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  smin = mesh->rphysical[eInterp*mesh->Np];
	  smax = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (snew - smin)/(smax - smin);
	  perp_index = n%mesh->Nq;
	  
	  rsdir = sdir;

	  break;
	case 33:
	  rold = mesh->rphysical[eAdj*mesh->Np + n];
	  sold = mesh->sphysical[eAdj*mesh->Np + n] - offset;

	  snew = -1*(sold + M_PI/2.);
	  rnew = -1*atan(tan(rold)/tan(sold));

	  rmin = mesh->rphysical[eOverlap*mesh->Np];
	  rmax = mesh->rphysical[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rmin - rnew > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (rnew - rmax > tol) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;

	  rmin = mesh->rphysical[eInterp*mesh->Np];
	  rmax = mesh->rphysical[eInterp*mesh->Np + mesh->Np - 1];
	  
	  par_loc = (rnew - rmin)/(rmax - rmin);
	  perp_index = n/mesh->Nq;
	  
	  rsdir = rdir;

	  break;	  
	}
      
      /*      if (rsdir == sdir && (rnew > M_PI/4 + tol || rnew < -M_PI/4 - tol || (rnew < M_PI/4 - offset - tol && rnew > -M_PI/4 + offset + tol))) printf("bad bounding %lf %d\n",rnew,faceHash);
      if (rsdir == rdir && (snew > M_PI/4 + tol || snew < -M_PI/4 - tol || (snew < M_PI/4 - offset - tol && snew > -M_PI/4 + offset + tol))) printf("bad bounding %lf %d\n",snew,faceHash);
      if (mesh->cubeFaceNumber[eInterp] != mesh->cubeFaceNumber[eOverlap]) printf("indexing error %d\n",faceHash);
      if (mesh->cubeDistance[eInterp] != 0) printf("interior element detected %d\n",faceHash);
      if (rsdir == sdir && smin >= smax) printf("element orientation mistake %d\n",faceHash);
      if (rsdir == rdir && rmin >= rmax) printf("element orientation mistake %d\n",faceHash);
      if (rsdir == sdir && (((sloc - 1) > tol) || (sloc < -tol))) {
	printf("placement error with sloc %lf, sabs %lf, rold %lf, sold %lf bounds %lf %lf,element %d and hash %d\n",sloc,snew,rold,sold,smin,smax,mesh->cubeFaceNumber[eOverlap],faceHash);
      }
      if (rsdir == rdir && (((rloc - 1) > tol) || (rloc < -tol))) {
	printf("placement error with rloc %lf, rabs %lf, rold %lf, sold %lf bounds %lf %lf,element %d and hash %d\n",rloc,rnew,rold,sold,rmin,rmax,eOverlap,faceHash);
	}*/
      
      mesh->rphysical[e*mesh->Np + n] = rnew;
      mesh->sphysical[e*mesh->Np + n] = snew;

      mesh->par_loc[i*mesh->Np + n] = par_loc;
      mesh->perp_index[i*mesh->Np + n] = perp_index;

      mesh->overlapDirection[i] = rsdir;
      mesh->eInterp[i*mesh->Np + n] = eInterp;
    }
  }
}
