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

      dfloat xsph, ysph, zsph, norm;

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
	break;
      case 2: //positive constant y
	xsph = sin(faceScale*xlin)*cos(faceScale*zlin);
	ysph = cos(faceScale*xlin)*cos(faceScale*zlin);
	zsph = sin(faceScale*zlin)*cos(faceScale*xlin);
	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);
	xsph /= norm;
	ysph /= norm;
	zsph /= norm;
	break;
      case 3: //negative constant x
	xsph = -1*cos(faceScale*ylin)*cos(faceScale*zlin);
	ysph = sin(faceScale*ylin)*cos(faceScale*zlin);
	zsph = sin(faceScale*zlin)*cos(faceScale*ylin);
	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);
	xsph /= norm;
	ysph /= norm;
	zsph /= norm;
	break;
      case 4: //negative constant y
	xsph = sin(faceScale*xlin)*cos(faceScale*zlin);
	ysph = -1*cos(faceScale*xlin)*cos(faceScale*zlin);
	zsph = sin(faceScale*zlin)*cos(faceScale*xlin);
	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);
	xsph /= norm;
	ysph /= norm;
	zsph /= norm;
	break;
      case 5: //positive constant z
	xsph = sin(faceScale*xlin)*cos(faceScale*ylin);
	ysph = sin(faceScale*ylin)*cos(faceScale*xlin);
	zsph = cos(faceScale*xlin)*cos(faceScale*ylin);
	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);
	xsph /= norm;
	ysph /= norm;
	zsph /= norm;
	break;
      case 6: //negative constant z
	xsph =  sin(faceScale*xlin)*cos(faceScale*ylin);
	ysph =  sin(faceScale*ylin)*cos(faceScale*xlin);
	zsph = -1*cos(faceScale*xlin)*cos(faceScale*ylin);
	norm = sqrt(xsph*xsph+ysph*ysph+zsph*zsph);
	xsph /= norm;
	ysph /= norm;
	zsph /= norm;
	break;
      }

      if (fabs(xsph*xsph + ysph*ysph + zsph*zsph - 1) > 1e-12) printf("error %d\n",mesh->cubeFaceNumber[e]);
	  
      // project to sphere
      mesh->x[cnt] = xsph; 
      mesh->y[cnt] = ysph; 
      mesh->z[cnt] = zsph;
      
      ++cnt;
    }
  }
}
