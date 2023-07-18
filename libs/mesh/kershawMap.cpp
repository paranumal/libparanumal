
#include "mesh.hpp"

namespace libp {

  // 1D transformation at the right boundary.
  dfloat right(const dfloat eps, const dfloat x){
    return (x <= 0.5) ? (2-eps) * x : 1 + eps*(x-1);
  }
  
  // 1D transformation at the left boundary
  dfloat left(const dfloat eps, const dfloat x){
    return 1-right(eps,1-x);
  }
  
  // Transition from a value of "a" for x=0, to a value of "b" for x=1.  Optionally
  // smooth -- see the commented versions at the end.
  dfloat step(int mapType, const dfloat a, const dfloat b, dfloat x)
  {
    if (x <= 0) return a;
    if (x >= 1) return b;

    switch(mapType){
    case 1:
      return a + (b-a) * (x);
    case 2:
      return a + (b-a) * (x*x*(3-2*x));
    default:
      return a + (b-a) * (x*x*x*(x*(6*x-15)+10));
    }

    return 0;
  }

  // 3D version of a generalized Kershaw mesh transformation, see D. Kershaw,
  // "Differencing of the diffusion equation in Lagrangian hydrodynamic codes",
  // JCP, 39:375–395, 1981.
  //
  // The input mesh should be Cartesian nx x ny x nz with nx divisible by 6 and
  // ny, nz divisible by 2.
  //
  // The eps parameters are in (0, 1]. Uniform mesh is recovered for epsy=1
  void coordMap2D(int mapType, const dfloat epsy, 
		  const dfloat x, const dfloat y, 
		  dfloat *X, dfloat *Y)
  {
    *X = x;

    int layer = x*6.0;
    dfloat lambda = (x-layer/6.0)*6;

    // The x-range is split in 6 layers going from left-to-left, left-to-right,
    // right-to-left (2 layers), left-to-right and right-to-right yz-faces.
    switch (layer)
      {
      case 0:
	*Y = left(epsy, y);
	break;
      case 1:
      case 4:
	*Y = step(mapType, left(epsy, y), right(epsy, y), lambda);
	break;
      case 2:
	*Y = step(mapType, right(epsy, y), left(epsy, y), lambda/2);
	break;
      case 3:
	*Y = step(mapType, right(epsy, y), left(epsy, y), (1+lambda)/2);
	break;
      default:
	*Y = right(epsy, y);
	break;
      }
  }

  // 3D version of a generalized Kershaw mesh transformation, see D. Kershaw,
  // "Differencing of the diffusion equation in Lagrangian hydrodynamic codes",
  // JCP, 39:375–395, 1981.
  //
  // The input mesh should be Cartesian nx x ny x nz with nx divisible by 6 and
  // ny, nz divisible by 2.
  //
  // The eps parameters are in (0, 1]. Uniform mesh is recovered for epsy=epsz=1.
  void coordMap3D(int mapType, const dfloat epsy, const dfloat epsz,
		  const dfloat x, const dfloat y, const dfloat z,
		  dfloat *X, dfloat *Y, dfloat *Z)
  {
    *X = x;

    int layer = x*6.0;
    dfloat lambda = (x-layer/6.0)*6;

    // The x-range is split in 6 layers going from left-to-left, left-to-right,
    // right-to-left (2 layers), left-to-right and right-to-right yz-faces.
    switch (layer)
      {
      case 0:
	*Y = left(epsy, y);
	*Z = left(epsz, z);
	break;
      case 1:
      case 4:
	*Y = step(mapType, left(epsy, y), right(epsy, y), lambda);
	*Z = step(mapType, left(epsz, z), right(epsz, z), lambda);
	break;
      case 2:
	*Y = step(mapType, right(epsy, y), left(epsy, y), lambda/2);
	*Z = step(mapType, right(epsz, z), left(epsz, z), lambda/2);
	break;
      case 3:
	*Y = step(mapType, right(epsy, y), left(epsy, y), (1+lambda)/2);
	*Z = step(mapType, right(epsz, z), left(epsz, z), (1+lambda)/2);
	break;
      default:
	*Y = right(epsy, y);
	*Z = right(epsz, z);
	break;
      }
  }
  

  
  void kershawMap2D(int mapType,  dfloat param1, dfloat &x, dfloat &y){

    dfloat xn = x + .5;
    dfloat yn = y + .5;
  
    dfloat Xn, Yn;
  
    // assume bi-unit box
    coordMap2D(mapType, param1, xn, yn, &Xn, &Yn);
  
    x = Xn-.5;
    y = Yn-.5;
  }
  
  
  void kershawMap3D(int mapType,  dfloat param1, dfloat &x, dfloat &y, dfloat &z){
    
    dfloat xn = x+.5;
    dfloat yn = y+.5;
    dfloat zn = z+.5;
  
    dfloat Xn, Yn, Zn;
  
    // assume bi-unit box
    coordMap3D(mapType, param1, param1, xn, yn, zn, &Xn, &Yn, &Zn);
  
    x = Xn-.5;
    y = Yn-.5;
    z = Zn-.5;
  }
}// end namespace content
