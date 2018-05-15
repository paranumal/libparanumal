#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

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

  mesh->r = (dfloat*) calloc(mesh->NgridElements*mesh->Np,sizeof(dfloat));
  mesh->s = (dfloat*) calloc(mesh->NgridElements*mesh->Np,sizeof(dfloat));

  mesh->rlocal = (dfloat*) calloc(mesh->NgridElements*mesh->Np,sizeof(dfloat));
  mesh->slocal = (dfloat*) calloc(mesh->NgridElements*mesh->Np,sizeof(dfloat));
  
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

      mesh->rlocal[cnt] = n%mesh->Nq;
      mesh->slocal[cnt] = n/mesh->Nq;

      mesh->r[cnt] = rlin;
      mesh->s[cnt] = slin;
      
      ++cnt;
    }
  }
}

void meshEquiSphericalExtensionQuad3D(mesh_t *mesh) {

  mesh->eInterp = (iint *) calloc((mesh->NgridElements - mesh->Nelements);
  
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
      mesh->rlocal[e*mesh->Np + n] = n%mesh->Nq;
      mesh->slocal[e*mesh->Np + n] = n/mesh->Nq;

      offset = M_PI/(2*mesh->edgeLength);
      
      iint faceHash = 6*face1 + face2;

      dfloat rold, sold, rnew, snew, rloc, sloc, rmin, smin, rmax, smax;
      iint rsdir;
      iint rdir = 0;
      iint sdir = 1;
      
      switch (faceHash)
	{
	case 6:
	case 13:
	case 20:
	case 3:
	  rold = mesh->r[eAdj*mesh->Np + n] + offset;
	  sold = mesh->s[eAdj*mesh->Np + n];

	  rnew = rold - M_PI/2;
	  snew = atan(tan(sold)/tan(rold));

	  smin = mesh->s[eOverlap*mesh->Np];
	  smax = mesh->s[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (snew < smin) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew > smax) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;

	  rloc = mesh->rlocal[eInterp*mesh->Np + n];
	  sloc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;

	  break;
	case 1:
	case 8:
	case 15:
	case 18:
	  rold = mesh->r[eAdj*mesh->Np + n] - offset;
	  sold = mesh->s[eAdj*mesh->Np + n];

	  rnew = rold + M_PI/2;
	  snew = atan(tan(sold)/tan(rold));

	  smin = mesh->s[eOverlap*mesh->Np];
	  smax = mesh->s[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (snew < smin) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 0];
	  else if (snew > smax) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 2];
	  else eInterp = eOverlap;

	  rloc = mesh->rlocal[eInterp*mesh->Np + n];
	  sloc = (snew - smin)/(smax - smin);
	  
	  rsdir = sdir;

	  break;
	case 4:
	case 30:
	  rold = mesh->r[eAdj*mesh->Np + n];
	  sold = mesh->s[eAdj*mesh->Np + n] - offset;

	  rnew = atan(tan(rold)/tan(sold));
	  snew = sold + M_PI/2;

	  rmin = mesh->r[eOverlap*mesh->Np];
	  rmax = mesh->r[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rnew < rmin) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (rnew > rmax) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  rloc = (rnew - rmin)/(rmax - rmin);
	  sloc = mesh->slocal[eInterp*mesh->Np + n];
	  
	  rsdir = rdir;

	  break;	  
	case 5:
	case 24:
	  rold = mesh->r[eAdj*mesh->Np + n];
	  sold = mesh->s[eAdj*mesh->Np + n] + offset;

	  rnew = atan(tan(rold)/tan(sold));
	  snew = sold - M_PI/2;

	  rmin = mesh->r[eOverlap*mesh->Np];
	  rmax = mesh->r[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rnew < rmin) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (rnew > rmax) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  rloc = (rnew - rmin)/(rmax - rmin);
	  sloc = mesh->slocal[eInterp*mesh->Np + n];
	  
	  rsdir = rdir;

	  break;
	case 16:
	case 32:
	  rold = mesh->r[eAdj*mesh->Np + n];
	  sold = mesh->s[(eAdj+1)*mesh->Np - n - 1] + offset;

	  rnew = atan(tan(rold)/tan(sold));
	  snew = sold - offset;

	  rmin = mesh->s[eOverlap*mesh->Np];
	  rmax = mesh->s[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rnew < rmin) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (rnew > rmax) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  rloc = (rnew - rmin)/(rmax - rmin);
	  sloc = mesh->slocal[(eInterp+1)*mesh->Np - n - 1];
	  
	  rsdir = rdir;

	  break;
	case 17:
	case 26:
	  rold = mesh->r[eAdj*mesh->Np + n];
	  sold = mesh->s[(eAdj+1)*mesh->Np - n - 1] - offset;

	  rnew = atan(tan(rold)/tan(sold));
	  snew = sold + M_PI/2;

	  rmin = mesh->s[eOverlap*mesh->Np];
	  rmax = mesh->s[eOverlap*mesh->Np + mesh->Np - 1];
	  
	  if (rnew < rmin) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 3];
	  else if (rnew > rmax) eInterp = mesh->EToE[eOverlap*mesh->Nfaces + 1];
	  else eInterp = eOverlap;

	  rloc = (rnew - rmin)/(rmax - rmin);
	  sloc = mesh->slocal[(eInterp+1)*mesh->Np - n - 1];
	  
	  rsdir = rdir;

	  break;
	case 25:

	case 10:

	case 22:

	case 27:

	case 11:

	case 31:

	case 23:

	case 33:
	  
	}
    }
  }
}
