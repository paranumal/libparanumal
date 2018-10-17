#include "cuda.h"
// risky
#define dfloat double 
#include "simpleRayTracer.h"

#define p_eps 1e-6

#define p_Nsamples 1

// ratio of importance in sampling primary ray versus random rays
#define p_primaryWeight 2.f

#define p_intersectDelta 0.1f

#define p_shadowDelta 0.15f
#define p_projectDelta 1e-2

#define p_maxLevel 5
#define p_maxNrays (2<<p_maxLevel)
#define p_apertureRadius 20.f
#define NRANDOM 10000

cudaEvent_t startTimer, endTimer;

void initTimer(){

  cudaEventCreate(&startTimer);
  cudaEventCreate(&endTimer);
}

void ticTimer(){

  cudaEventRecord(startTimer);
}

void tocTimer(const char *message){

  cudaEventRecord(endTimer);
  cudaEventSynchronize(endTimer);

  float elapsed;
  cudaEventElapsedTime(&elapsed, startTimer, endTimer);

  printf("Kernel %s took %g seconds\n", message, elapsed/1000.);
}


__device__ bbox_t createBoundingBoxSphere(sphere_t &sphere);

__host__ __device__ dfloat clamp(dfloat x, dfloat xmin, dfloat xmax){

  x = min(x, xmax);
  x = max(x, xmin);
  
  return x;
}


__host__ __device__ int iclamp(dfloat x, dfloat xmin, dfloat xmax){

  x = min(x, xmax);
  x = max(x, xmin);
  
  return floor(x);
}



__forceinline__ __host__ __device__  vector_t sensorLocation(const int NI,
							     const int NJ,
							     const int I,
							     const int J,
							     const sensor_t &sensor){
  
  vector_t sensorX = sensor.eyeX;
  
  dfloat r = I/(dfloat)(NI-1);
  dfloat s = J/(dfloat)(NJ-1);
  
  r = (r-0.5f)*sensor.Ilength;
  s = (s-0.5f)*sensor.Jlength;

  vector_t sensorNormal =
    vectorCrossProduct(sensor.Idir, sensor.Jdir);
  
  sensorX.x += sensorNormal.x*sensor.offset;
  sensorX.y += sensorNormal.y*sensor.offset;
  sensorX.z += sensorNormal.z*sensor.offset;

  sensorX.x += r*sensor.Idir.x;
  sensorX.y += r*sensor.Idir.y;
  sensorX.z += r*sensor.Idir.z;
  
  sensorX.x += s*sensor.Jdir.x;
  sensorX.y += s*sensor.Jdir.y;
  sensorX.z += s*sensor.Jdir.z;

  return sensorX;
}

__host__ __device__  void sensorMultipleLocations(const int NI,
						  const int NJ,
						  const int I,
						  const int J,
						  const sensor_t &sensor, 
						  vector_t *sensorsX){
  
  for(int samp=0;samp<p_Nsamples;++samp){
    sensorsX[samp] = sensor.eyeX;
    
    dfloat r = I/(dfloat)(NI-1);
    dfloat s = J/(dfloat)(NJ-1);

    dfloat theta = 2.f*M_PI*samp/(dfloat)p_Nsamples;
    
    // circle of samples around sensor pixel
    dfloat delta = .5; // scatter pixel radius
    r = (r-0.5f+delta*cosf(theta)/NI)*sensor.Ilength;
    s = (s-0.5f+delta*sinf(theta)/NJ)*sensor.Jlength;

    vector_t sensorNormal = vectorCrossProduct(sensor.Idir, sensor.Jdir);

    sensorsX[samp].x += sensorNormal.x*sensor.offset;
    sensorsX[samp].y += sensorNormal.y*sensor.offset;
    sensorsX[samp].z += sensorNormal.z*sensor.offset;
    
    sensorsX[samp].x += r*sensor.Idir.x;
    sensorsX[samp].y += r*sensor.Idir.y;
    sensorsX[samp].z += r*sensor.Idir.z;

    sensorsX[samp].x += s*sensor.Jdir.x;
    sensorsX[samp].y += s*sensor.Jdir.y;
    sensorsX[samp].z += s*sensor.Jdir.z;
  }

}


__host__ __device__  vector_t vectorCreate(dfloat x, dfloat y, dfloat z){
  vector_t v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

/* Subtract two vectors and return the resulting vector_t */
__host__ __device__  vector_t vectorSub(const vector_t v1, const vector_t v2){
  return vectorCreate(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

/* Multiply two vectors and return the resulting scalar (dot product) */
__host__ __device__  dfloat vectorDot(const vector_t v1, const vector_t v2){
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

/* Entrywise multiply two vectors and return the resulting vector */
__host__ __device__  vector_t vectorDotMultiply(const vector_t v1, const vector_t v2){
  return vectorCreate(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

/* Entrywise divison of  two vectors and return the resulting vector */
__host__ __device__  vector_t vectorDotDivide(const vector_t v1, const vector_t v2){
  return vectorCreate(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}

__host__ __device__  vector_t vectorCrossProduct(const vector_t v1, const vector_t v2){
  return vectorCreate(v1.y*v2.z-v1.z*v2.y,
		      v1.z*v2.x-v1.x*v2.z,
		      v1.x*v2.y-v1.y*v2.x);
}


/* Calculate Vector_T x Scalar and return resulting Vector*/ 
__host__ __device__  vector_t vectorScale(const dfloat c, const vector_t v){
  return vectorCreate(v.x * c, v.y * c, v.z * c);
}

/* Add two vectors and return the resulting vector_t */
__host__ __device__  vector_t vectorAdd(const vector_t v1, const vector_t v2){
  return vectorCreate(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

__host__ __device__  dfloat vectorTripleProduct(const vector_t a, const vector_t b, const vector_t c){

  const vector_t aXb = vectorCrossProduct(a, b);
  
  return vectorDot(aXb, c); 
}

// assume b is unit vector
__host__ __device__ vector_t vectorOrthogonalize(const vector_t a, const vector_t b){

  dfloat adotb = vectorDot(a, b);

  return  vectorSub(a, vectorScale(adotb, b));
}

__host__ __device__ dfloat vectorNorm(const vector_t a){
  return  sqrt(vectorDot(a,a));
}

// return orthonormalized vector
__host__ __device__ vector_t vectorNormalize(const vector_t a){

  dfloat d = vectorNorm(a);
  if(d)
    return vectorScale(1./d, a);
  else
    return vectorCreate(0,0,0);
}

// https://www.scratchapixel.com/code.php?id=10&origin=/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes
// roots of a*t^2 + 2*b*t + c = 0
__forceinline__  __host__ __device__ bool solveQuadratic(const dfloat &a, const dfloat &b, const dfloat &c, dfloat &x0, dfloat &x1){
  
  dfloat discr = b * b - a * c;

  if (discr < 0) return false;
  else if (discr == 0) {
    x0 = x1 = - b / a;
  }
  else {
    dfloat sqrtdiscr = sqrt(discr);
    dfloat q = (b > 0) ?
      -(b + sqrtdiscr) :
      -(b - sqrtdiscr);
    x0 = q / a;
    x1 = c / q;
  }

  dfloat xmin = min(x0, x1);
  dfloat xmax = max(x0, x1);
  x0 = xmin;
  x1 = xmax;
  
  
  return true; 
}

/* Check if the ray and triangle intersect */
__forceinline__  __host__ __device__ bool intersectRayTriangle(const ray_t &r, const triangle_t &tri, dfloat *t){

  // TW: unused fudge factor
  dfloat delta = 0; 
  
  bool retval = false;
 
  vector_t B1 = vectorSub(tri.vertices[2], tri.vertices[0]);
  vector_t B2 = vectorSub(tri.vertices[2], tri.vertices[1]);
  vector_t B3 = r.dir;

   vector_t R = vectorSub(tri.vertices[2], r.start);

  dfloat J = vectorTripleProduct(B2, B3, B1);
  
  dfloat L1 = vectorTripleProduct(B2, B3, R);
  if(L1<delta*J) return false;

  dfloat L2 = vectorTripleProduct(B3, B1, R);
  if(L2<delta*J || L1+L2>J*(1+delta)) return false;

  dfloat t0 = vectorTripleProduct(B1, B2, R)/J;

  /* Verify t1 larger than 0 and less than the original t */
  // TW: FUDGE FACTOR
  if((t0 > p_intersectDelta) && (t0 < *t)){
    *t = t0;
    retval = true;
  }

  return retval;
}

/* Check if the ray and triangle intersect */
__forceinline__  __host__ __device__ bool intersectRayRectangle(const ray_t &r, const rectangle_t &rect, dfloat *t){

  vector_t C  = rect.center;
  vector_t A1 = rect.axis[0];
  vector_t A2 = rect.axis[1];
  dfloat   L1 = rect.length[0];
  dfloat   L2 = rect.length[1];

  // n = A1 x A2
  // (s + t*d - C).n  = 0
  // t = (C - s).n/(d.n)

  vector_t n = vectorCrossProduct(A1, A2);
  dfloat  t0 = vectorDot(vectorSub(C,r.start), n)/vectorDot(r.dir, n);

  // intersection behind start of ray
  if(t0<0 || t0>*t) return false;

  // X = s + t*d - C
  vector_t X = vectorAdd(vectorSub(r.start,C), vectorScale(t0, r.dir));
  
  dfloat h1 = vectorDot(A1, X)+0.5*L1; // shift
  if(h1<0 || h1>L1) return false;
  
  dfloat h2 = vectorDot(A2, X)+0.5*L2; // shift
  if(h2<0 || h2>L1) return false;

  // success
  *t = t0;
  
  return true;
}



/* Check if the ray and sphere intersect */
__forceinline__ __host__ __device__ bool intersectRaySphere(const ray_t &r, const sphere_t &s, dfloat *t){
	
  bool retval = false;

  /* A = d.d, the vector_t dot product of the direction */
  dfloat A = vectorDot(r.dir, r.dir); 
	
  /* We need a vector_t representing the distance between the start of 
   * the ray and the position of the circle.
   * This is the term (p0 - c) 
   */
  vector_t dist = vectorSub(r.start, s.pos);
	
  /* 2d.(p0 - c) */  
  dfloat B = 2.f * vectorDot(r.dir, dist);
	
  /* (p0 - c).(p0 - c) - r^2 */
  dfloat C = vectorDot(dist, dist) - (s.radius * s.radius);

  /* find roots of quadratic */
  dfloat t0, t1;
  if(solveQuadratic(A,0.5*B,C,t0,t1)){
    if((t0 > p_intersectDelta) && (t0 < *t)){
      *t = t0;
      retval = true;
    }else
      retval = false;
  }else{
    retval = false;
  }

  return retval;
}

/* Check if the ray and sphere intersect */
__forceinline__ __host__ __device__ bool intersectRayEllipsoid(const ray_t &r, const ellipsoid_t &sh, dfloat *t){
	
  bool retval = false;

  /* 
     R = vector radii
     ((p.x-c.x)/R.x)^2 + ((p.y-c.y)/R.y)^2 + ((p.z-c.z)/R.z)^2 = 1;

     p.x = s.x + t*d.x

     s~ = s-c

     ((s~.x+t*d.x)/R.x)^2 + ((s~.y+t*d.y)/R.y)^2 + ((s~.z+t*d.z)/R.z)^2 = 1;

     t^2 *   ( (d.x/R.x)^2 + (d.y/R.y)^2 + (d.z/R.z)^2 
    +t   *2* ( (d.x*s~.x/R.x^2) + (d.y*s~.y/R.y^2) + (d.z*s~.z/R.z^2) )
    +        ( (s~.x/R.x)^2 + (s~.y/R.y)^2 + (s~.z/R.z)^2 - 1 ) = 0
  */

  vector_t s = r.start;
  vector_t d = r.dir;
  vector_t c = sh.pos;
  vector_t invR = sh.invRadii; 
  vector_t st = vectorSub(s, c);
  
  /* A = d.d, the vector_t dot product of the direction */
  vector_t dIR  = vectorDotMultiply(d,  invR);
  vector_t stIR = vectorDotMultiply(st, invR);

  dfloat A = vectorDot(dIR, dIR);
  dfloat B = vectorDot(dIR, stIR);
  dfloat C = vectorDot(stIR, stIR) - 1.f;  

  /* find roots of quadratic */
  dfloat t0, t1;
  if(solveQuadratic(A,B,C,t0,t1)){
    if((t0 > p_intersectDelta) && (t0 < *t)){
      *t = t0;
      retval = true;
    }else
      retval = false;
  }else{
    retval = false;
  }

  return retval;
}



__forceinline__ __host__ __device__ bool intersectRayCone(const ray_t &r, const cone_t &cone, dfloat *t){

  bool retval = false;

  /* 

     cone-ray intersection tests:
     
     | pos + t*dir - (vertex + axis*h) | = R*h/H

  */

  vector_t p = r.start;
  vector_t d = r.dir;

  vector_t v = cone.vertex;
  vector_t a = cone.axis;
  dfloat   R = cone.radius;
  dfloat   H = cone.height;
  
  dfloat alpha = (R/H)*(R/H);

  // p + t*d - (v + h*a) orth a
  // h = (p-v + t*d).a
  // if h>=0
  // | p + t*d - (v + h*a)| = alpha*(p-v + t*d).a 
  // |(p-v) - ((p-v).a)*a + t*d - t*(a.d)*a | = alpha*(p-v+t*d).a
  // | pminusvPerp + t*dPerp| = alpha*( (p-v).a + t*d.a)
  // 
  
  dfloat adotd = vectorDot(a,d);
  vector_t dPerp = vectorSub(d, vectorScale(adotd, a));

  vector_t pminusv = vectorSub(p, v);
  dfloat tmp = vectorDot(a,pminusv);
  
  vector_t pminusvPerp = vectorSub(pminusv, vectorScale(tmp, a));
  
  dfloat A = vectorDot(dPerp, dPerp)             - alpha*adotd*adotd;
  dfloat B = vectorDot(dPerp, pminusvPerp)       - alpha*adotd*tmp;
  dfloat C = vectorDot(pminusvPerp, pminusvPerp) - alpha*tmp*tmp;

  /* find roots of quadratic */
  dfloat t0, t1;
  if(solveQuadratic(A,B,C,t0,t1)){

    // cone is behind ray
    if(t0<0 && t1<0)
      return false;

    // check location along axis
    const dfloat h0 = tmp + t0*adotd;
    const dfloat h1 = tmp + t1*adotd;

    const int valid0 = ((h0>0) && (h0<H));
    const int valid1 = ((h1>0) && (h1<H));

    if(!valid0 && !valid1) // out of range
      return false;
    else if(valid0 && valid1){ // both viable
      if(t0 > t1){ // nearest
	t0 = t1;
      }
    }
    else if(valid1){
      t0 = t1;
    }
    
    if((t0 > p_intersectDelta) && (t0 < *t)){
      *t = t0;
      retval = true;
    }else
      retval = false;
  }else{
    retval = false;
  }
  
  return retval;
}

__forceinline__ __host__ __device__ bool intersectRayDisk(const ray_t &r, const disk_t &disk, dfloat *t){

  vector_t s = r.start;
  vector_t d = r.dir;
  vector_t n = disk.normal;
  vector_t c = disk.center;
  
  // intersection with plane
  dfloat ndotd = vectorDot(n, d);

  // (s + t*d -c ).n = 0
  dfloat  t0 = vectorDot(vectorSub(c,s), n)/ndotd;

  // intersection behind start of ray
  if(t0<0 || t0>*t) return false;

  vector_t p = vectorAdd(s, vectorScale(t0, d));
  vector_t v = vectorSub(p, c);

  dfloat R2 = vectorDot(v,v);

  if(R2>=(disk.radius*disk.radius)-p_intersectDelta)
    return false;

  if(t0>*t) return false;
  
  *t = t0;
  return true;
}

__forceinline__ __host__ __device__ bool intersectRayCylinder(const ray_t &r, const cylinder_t &cylinder, dfloat *t){

  bool retval = false;

  /* 

     cylinder-ray intersection tests:
     
     | p  + t*d - (c+h*a) | = R
     h = (p+t*d-c).a
     0<= h <=H

  */

  vector_t p = r.start;
  vector_t d = r.dir;

  vector_t c = cylinder.center;
  vector_t a = cylinder.axis;
  dfloat   R = cylinder.radius;
  dfloat   H = cylinder.height;
  dfloat adotd = vectorDot(a,d);
  vector_t dPerp = vectorSub(d, vectorScale(adotd, a));

  vector_t pminusc = vectorSub(p, c);
  dfloat tmp = vectorDot(a,pminusc);
  vector_t pminuscPerp = vectorSub(pminusc, vectorScale(tmp, a));
  
  dfloat A = vectorDot(dPerp, dPerp);
  dfloat B = vectorDot(dPerp, pminuscPerp);
  dfloat C = vectorDot(pminuscPerp, pminuscPerp) - R*R;

#if 1
  // prone to acne (FP32)
  dfloat t0, t1;
  if(solveQuadratic(A,B,C,t0,t1)){

    // cylinder is behind ray
    if(t0<=0 && t1<=0)
      return false;
    
    dfloat h0 = tmp + t0*adotd;
    dfloat h1 = tmp + t1*adotd;

    int valid0 = ((h0>0) && (h0<H));
    int valid1 = ((h1>0) && (h1<H));

    if(!valid0 && !valid1){      
      return false;
    }
    else if(valid0 && valid1){
      if(t0 > t1){
	t0= t1;
      }
    }
    else if(valid1){
      t0 = t1;
    }

    // TW: FUDGE FACTOR (was 1e-3)
    if((t0 > p_intersectDelta) && (t0< ((*t)))){// weakened this test
      *t = t0;
      retval = true;
    }else
      retval = false;
  }else{
    retval = false;
  }
#else
  // prone to acne (FP32)
  dfloat discr = B*B-A*C;

  // TW: UNUSED FUDGE FACTOR
  dfloat delta = p_intersectDelta; // need large tolerance
  if(discr<=delta)
    retval = false;
  else{
    dfloat sqrtdiscr = sqrtf(discr);
    dfloat A2 = A*A;
    dfloat t0A2 = (-B + sqrtdiscr)*A;
    dfloat t1A2 = (-B - sqrtdiscr)*A;

    if(t0A2<=delta*A2 && t1A2<=delta*A2) return false;
    
    dfloat h0A2 = tmp*A2 + t0A2*adotd;
    dfloat h1A2 = tmp*A2 + t1A2*adotd;

    int valid0 = ((h0A2>delta*A2) && (h0A2<H*A2-delta*A2));
    int valid1 = ((h1A2>delta*A2) && (h1A2<H*A2-delta*A2));

    if(!valid0 && !valid1)
      return false;
    else if(valid0 && valid1){
      if(t0A2 > t1A2){
	t0A2 = t1A2;
      }
    }
    else if(valid1){
      t0A2 = t1A2;
    }

    // TW: FUDGE FACTOR (was 1e-3)
    if((t0A2 > p_intersectDelta*A2) && (t0A2 < ((*t)*p_intersectDelta))){// weakened this test
      *t = t0A2/A2;
      retval = true;
    }else
      retval = false;
  }
#endif
  
  return retval;
}

__host__ __device__ bool intersectPointGridCell(const grid_t &grid,
						const vector_t p,
						const int cellI,
						const int cellJ,
						const int cellK){

  
  
  if(p.x<=grid.xmin+(cellI  )*grid.dx) return false;
  if(p.x> grid.xmin+(cellI+1)*grid.dx) return false;

  if(p.y<=grid.ymin+(cellJ  )*grid.dy) return false;
  if(p.y> grid.ymin+(cellJ+1)*grid.dy) return false;
  
  if(p.z<=grid.zmin+(cellK  )*grid.dz) return false;
  if(p.z> grid.zmin+(cellK+1)*grid.dz) return false;  
  
  return true;
}

__host__ __device__ bool intersectRayBox(ray_t &r, const bbox_t &bbox, unsigned int &face){

  vector_t d = r.dir;
  vector_t s = r.start;
  vector_t invd = r.invDir;
  
  dfloat mint = 20000;
  face = 0;
  
  if(d.x>0){ // face 2
    dfloat newt = (bbox.xmax-s.x)*invd.x; // d.x > 0
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.x<0){ // face 4
    // s.x + newt*d.x = bbox.xmin
    dfloat newt = (bbox.xmin-s.x)*invd.x;
    if(newt>0){
      mint = min(mint, newt);
    }
  }
  
  if(d.y>0){ // face 3
    dfloat newt = (bbox.ymax-s.y)*invd.y;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.y<0){ // face 1
    dfloat newt = (bbox.ymin-s.y)*invd.y;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.z>0){ // face 5
    dfloat newt = (bbox.zmax-s.z)*invd.z;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.z<0){ // face 0
    dfloat newt = (bbox.zmin-s.z)*invd.z;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  face = 0;
  
  if(d.x>0){ // face 2
    dfloat newt = (bbox.xmax-s.x)*invd.x;
    if(newt>0 && newt<=mint)
      face |= 4;
  }

  if(d.x<0){ // face 4
    dfloat newt = (bbox.xmin-s.x)*invd.x;
    if(newt>0 && newt<=mint)
      face |= 16;
  }

  if(d.y>0){ // face 3
    dfloat newt = (bbox.ymax-s.y)*invd.y;
    if(newt>0 && newt<=mint)
      face |= 8;
  }

  if(d.y<0){ // face 1
    dfloat newt = (bbox.ymin-s.y)*invd.y;
    if(newt>0 && newt<=mint)
      face |= 2;
  }

  if(d.z>0){ // face 5
    dfloat newt = (bbox.zmax-s.z)*invd.z;
    if(newt>0 && newt<=mint)
      face |= 32;
  }

  if(d.z<0){ // face 0
    dfloat newt = (bbox.zmin-s.z)*invd.z;
    if(newt>0 && newt<=mint)
      face |= 1;
  }
  
  if(face>0){
    r.start = vectorAdd(s, vectorScale(mint, d));
    return true;
  }

  return false;
}

__forceinline__ __host__ __device__ bool intersectRayShape(const ray_t &r, const shape_t &s, dfloat *t){

  switch(s.type){
  case SPHERE:   return intersectRaySphere   (r, s.sphere,   t);
  case CONE:     return intersectRayCone     (r, s.cone,     t);
  case DISK:     return intersectRayDisk     (r, s.disk,     t);
  case CYLINDER: return intersectRayCylinder (r, s.cylinder, t);
  case IMAGE:
  case RECTANGLE:return intersectRayRectangle(r, s.rectangle, t);
  case TRIANGLE: return intersectRayTriangle (r, s.triangle,  t); 
  case ELLIPSOID:return intersectRayEllipsoid(r, s.ellipsoid, t);    
  }

  return false;
  
}


__forceinline__ __host__ __device__ vector_t computeNormal(const vector_t &v, const shape_t &s){

  vector_t n = vectorCreate(0,0,0);
  
  /* Find the normal for this new vector_t at the point of intersection */
  switch(s.type){
  case SPHERE:
    {
      n = vectorSub(v, s.sphere.pos);
      break;
    }
  case ELLIPSOID:
    {
      vector_t vMs = vectorSub(v, s.ellipsoid.pos);

      // f = (v-c).^2./(radii.^2) - 1 => n = grad f
      n = vectorDotMultiply(vMs, vectorDotMultiply(s.ellipsoid.invRadii, s.ellipsoid.invRadii));

      break;
    }
    
  case TRIANGLE:
    {
      vector_t a = vectorSub(s.triangle.vertices[2], s.triangle.vertices[0]); 
      vector_t b = vectorSub(s.triangle.vertices[1], s.triangle.vertices[0]);
      
      n = vectorCrossProduct(a, b);
      break;
    }
  case CONE:
    {
      // n = (v-vertex) x ( a x (v-vertex) )
      vector_t vMinusVertex = vectorSub(v, s.cone.vertex);

      // axis location
      dfloat H = s.cone.height;
      dfloat z = vectorDot(vMinusVertex, s.cone.axis);
      
      // problematic if axis is parallel to v-Vertex
      if(z>p_projectDelta && z<H-p_projectDelta)
	n = vectorCrossProduct( vMinusVertex, vectorCrossProduct(s.cone.axis, vMinusVertex));

      break;
    }
  case DISK:
    {
      vector_t vMc = vectorSub(v, s.disk.center);
      
      dfloat   R = s.disk.radius;

      vector_t tmp = vectorOrthogonalize(vMc, s.disk.normal);
      
      dfloat z = vectorNorm(tmp);

      if(z<R-p_projectDelta)
	n =  s.disk.normal;

      break;
    }
  case CYLINDER:
    {
      // z = (v - c).a  => clamp

      vector_t vMc = vectorSub(v, s.cylinder.center);

      dfloat   H = s.cylinder.height;
      dfloat z = vectorDot(vMc, s.cylinder.axis);

      if(z>p_projectDelta && z<H-p_projectDelta)
	n = vectorOrthogonalize(vMc, s.cylinder.axis);
      
      break;
    }
  case IMAGE:
  case RECTANGLE:
    {
#if 0
      vector_t C  = s.rectangle.center;
      vector_t A1 = s.rectangle.axis[0];
      vector_t A2 = s.rectangle.axis[1];
      dfloat   L1 = s.rectangle.length[0];
      dfloat   L2 = s.rectangle.length[1];
      
      
      // X = v - C
      vector_t X = vectorSub(v, C);
      
      dfloat h1 = vectorDot(A1, X)+0.5*L1; // shift
      dfloat h2 = vectorDot(A2, X)+0.5*L2; // shift
#endif
      
      n = vectorCrossProduct(s.rectangle.axis[0], s.rectangle.axis[1]);
      break;
    }
  }

  // normalize when normal is not degenerate
  dfloat tmp = vectorNorm(n);
  if(tmp)
    n = vectorScale(1./tmp, n);
  
  return n;
}

__forceinline__ __host__ __device__ material_t computeMaterial(const int Nmaterials, const material_t *materials,
							     const vector_t &v, const shape_t &s){
  material_t m;
  
  switch(s.type){
  case TRIANGLE:
    {

      // v = L1*v1 + L2*v2 + (1-L1-L2)*v3 + N1*n
      // [ v3-v1 v2-v1 -N1][L1;L2;n] = [v3-v]
      
      vector_t B1 = vectorSub(s.triangle.vertices[2], s.triangle.vertices[0]);
      vector_t B2 = vectorSub(s.triangle.vertices[1], s.triangle.vertices[0]);
      vector_t B3 = vectorCrossProduct(B1,B2);
      
      vector_t R = vectorSub(s.triangle.vertices[2], v);

      dfloat J  = vectorTripleProduct(B2, B3, B1);      
      dfloat L1 = vectorTripleProduct(B2, B3, R)/J;
      dfloat L2 = vectorTripleProduct(B3, B1, R)/J;
      dfloat N1 = vectorTripleProduct(B1, B2, R)/J;      

      dfloat Iq = L1*s.triangle.q[0] + L2*s.triangle.q[1] + (1-L1-L2)*s.triangle.q[2];

//      dfloat maxIq = 2.5, minIq =.05;
      dfloat maxIq = 3, minIq = -3;

      Iq = (Iq-minIq)/(maxIq-minIq);
#if 0
      if(Iq<0 || Iq>1){
	m.diffuse.red = 1;
	m.diffuse.green = 1;
	m.diffuse.blue = 1;
	m.reflection = 0;
	m.eta = 1.;
	m.refraction = 1;
	
	m.info.refractor = 1;
	m.info.reflector = 0;
	m.info.emitter = 0;
      }
      else{
      }
#endif	
      Iq = clamp(Iq, 0, 1);

#if 0
      dfloat redIq   = (Iq<1./3.) ? 3.*Iq:0; // reverse ?
      dfloat greenIq = (1./3<=Iq && Iq<2./3.) ? 3.*(Iq-1./3):0;
      dfloat blueIq  = (2./3<=Iq) ? 3.*(Iq-2./3):0;
#else
      dfloat redIq = 0, greenIq = 0, blueIq = 0;

      if(Iq<1/3.) redIq = 3*Iq;
      else if(Iq>=2./3) blueIq = 3.*(Iq-2./3);
      else{
	redIq = 1;
	greenIq = 1;
	blueIq = 1;
      }
	
      
#endif     
      m.diffuse.red   = redIq;
      m.diffuse.green = greenIq;
      m.diffuse.blue  = blueIq;
      m.reflection = 0.05;
      m.eta = 1.;
      m.refraction = 0.01;
      
      m.info.refractor = 0;
      m.info.reflector = 1;
      m.info.emitter = 0;
          
      break;
    }
  case SPHERE:
  case ELLIPSOID:
  case DISK:
  case CYLINDER:
  case CONE:
    {
      m = materials[s.material];
      break;
    }
#if 0
  case CYLINDER:
    {
      vector_t c = s.cylinder.center;
      vector_t a = s.cylinder.axis;
      vector_t vMc = vectorSub(v, c);
      dfloat H = s.cylinder.height;
      dfloat h = vectorDot(vMc, a);

      int i = (int) (8.f*(h/H)); // checkerboard material selector
      
      int idM = 10*((i%2)); // 1 if either i is odd or j is even
      m = materials[idM];
      break;
    }
  case CONE:
    {
      vector_t c = s.cone.vertex;
      vector_t a = s.cone.axis;
      vector_t vMc = vectorSub(v, c);
      dfloat H = s.cone.height;
      dfloat h = vectorDot(vMc, a);

      int i = (int) (8.f*(h/H)); // checkerboard material selector
      
      int idM = 20*((i%2)); // 1 if either i is odd or j is even
      m = materials[idM];
      break;
    }
#endif    
  case RECTANGLE:
    {
      if(s.material>=0)
	m = materials[s.material];
      else{
	vector_t C  = s.rectangle.center;
	vector_t A1 = s.rectangle.axis[0];
	vector_t A2 = s.rectangle.axis[1];
	dfloat   L1 = s.rectangle.length[0];
	dfloat   L2 = s.rectangle.length[1];
	
	
	// X = v - C
	vector_t X = vectorSub(v, C);
	
	dfloat h1 = vectorDot(A1, X)+0.5*L1; // shift
	dfloat h2 = vectorDot(A2, X)+0.5*L2; // shift
	
	int i = (int) (8.f*(h1/L1)); // checkerboard material selector
	int j = (int) (8.f*(h2/L2));
	
	int idM = ((i%2) ^ ((j+1)%2)); // 1 if either i is odd or j is even
	//      printf("i=%d, j=%d, h1=%g, h2=%g, L1=%g, L2=%g, idM = %d\n", i, j, h1, h2, L1, L2, idM);
	m = materials[idM];
      }
      break;
    }
  case IMAGE:
    {
      vector_t C  = s.rectangle.center;
      vector_t A1 = s.rectangle.axis[0];
      vector_t A2 = s.rectangle.axis[1];
      dfloat   L1 = s.rectangle.length[0];
      dfloat   L2 = s.rectangle.length[1];
      const unsigned char *img = s.rectangle.image;
      int      NI  = s.rectangle.NI;
      int      NJ  = s.rectangle.NJ;
      
      // X = v - C
      vector_t X = vectorSub(v, C);
      
      dfloat h1 = vectorDot(A1, X)+0.5*L1; // shift
      dfloat h2 = vectorDot(A2, X)+0.5*L2; // shift

      // nearest neighbor interpolation
      dfloat  i = iclamp(NI*h1/L1, 0, NI-1); 
      dfloat  j = iclamp(NJ*h2/L2, 0, NJ-1);
      int idM = (NI-1-i) + j*NI;

      m.diffuse.red   = img[idM*3 + 0]/256.f;
      m.diffuse.green = img[idM*3 + 1]/256.f;
      m.diffuse.blue  = img[idM*3 + 2]/256.f;

      m.reflection = 1;
      m.refraction = 1;
      m.info.refractor = 0;
      m.info.reflector = 0;
      m.info.emitter = 1;

      break;
    }
  }

  return m;
  
}


// grid search
__host__ __device__ bool gridRayIntersectionSearch(ray_t r,
						   const int Nshapes, const shape_t *shapes, const  grid_t &grid,
						   dfloat *t, int &currentShape){

  // is start of ray in a grid cell ?
  vector_t s = r.start; // will modify ray through s
  vector_t d = r.dir;
  vector_t invd;
  if(d.x)  invd.x = 1.f/d.x;
  if(d.y)  invd.y = 1.f/d.y;
  if(d.z)  invd.z = 1.f/d.z;
    
  // if ray is outside grid then project onto grid
  if(s.x<grid.xmin){
    if(d.x<=0) return false; // pointing away or grazing from grid
    dfloat t0 = -(s.x-grid.xmin)*invd.x;
    s.x = grid.xmin;
    s.y += t0*d.y;
    s.z += t0*d.z;
  }

  if(s.x>grid.xmax){
    if(d.x>=0) return false; 
    dfloat t0 = -(s.x-grid.xmax)*invd.x;
    s.x = grid.xmax;
    s.y += t0*d.y;
    s.z += t0*d.z;
  }

  if(s.y<grid.ymin){
    if(d.y<=0) return false; // pointing away or grazing from grid
    dfloat t0 = -(s.y-grid.ymin)*invd.y;
    s.y = grid.ymin;
    s.x += t0*d.x;
    s.z += t0*d.z;
  }

  if(s.y>grid.ymax){
    if(d.y>=0) return false; 
    dfloat t0 = -(s.y-grid.ymax)*invd.y;
    s.y = grid.ymax;
    s.x += t0*d.x;
    s.z += t0*d.z;
  }

  if(s.z<grid.zmin){
    if(d.z<=0) return false; // pointing away or grazing from grid
    dfloat t0 = -(s.z-grid.zmin)*invd.z;
    s.z = grid.zmin;
    s.x += t0*d.x;
    s.y += t0*d.y;
  }

  if(s.z>grid.zmax){
    if(d.z>=0) return false; 
    dfloat t0 = -(s.z-grid.zmax)*invd.z;
    s.z = grid.zmax;
    s.x += t0*d.x;
    s.y += t0*d.y;
  }
  
  // now the ray start must be on the surface of the grid or in a cell

  int cellI = iclamp((s.x-grid.xmin)*grid.invdx,0,grid.NI-1); // assumes grid.NI
  int cellJ = iclamp((s.y-grid.ymin)*grid.invdy,0,grid.NJ-1);
  int cellK = iclamp((s.z-grid.zmin)*grid.invdz,0,grid.NK-1);
  
  ray_t newr = r;
  newr.start = s;
  newr.invDir = invd;
  
  currentShape = -1;
  
  do{
    int cellID = cellI + grid.NI*cellJ + grid.NI*grid.NJ*cellK;
    
    *t = 20000; // TW ?

    int start = grid.c_boxStarts[cellID];
    int end   = grid.c_boxStarts[cellID+1];
    for(int offset=start;offset<end;++offset){
      const int obj = grid.c_boxContents[offset];
      const shape_t shape = shapes[obj];
      if(intersectRayShape(r, shape, t)){
	vector_t intersect = vectorAdd(r.start, vectorScale(*t, r.dir));
	
	if(intersectPointGridCell(grid, intersect, cellI, cellJ, cellK)){
	  currentShape = obj;
	}
      }
    }
    
    if(currentShape != -1){
      return true;
    }
    
    unsigned int face = 0;

    // find faces that ray passes through
    intersectRayBox(newr,grid.c_bboxes[cellID], face);
    
    if(face&1) --cellK; // face 0
    if(face&2) --cellJ; // face 1
    if(face&4) ++cellI; // face 2
    if(face&8) ++cellJ; // face 3
    if(face&16) --cellI;// face 4
    if(face&32) ++cellK;// face 5

    if(face==0){
      break;
    }
  }while(cellI>=0 && cellI<grid.NI &&
	 cellJ>=0 && cellJ<grid.NJ &&
	 cellK>=0 && cellK<grid.NK);

  return false;
}

__device__ colour_t trace(const grid_t grid,
			  const int Nshapes,
			  const shape_t *shapes,
			  const int Nlights,
			  const light_t *lights,
			  const int Nmaterials,
			  const material_t *materials,
			  ray_t  r,
			  int    level,
			  dfloat coef,
			  colour_t bg){
  
  colour_t black;
  black.red = 0;
  black.green = 0;
  black.blue = 0;

  // initialize color as black
  colour_t c = black;
  
  int Nrays = 0, rayID = 0;
  ray_t rayStack[p_maxNrays];

  // add initial ray to stack
  rayID = 0;
  r.level = 0;
  r.coef = coef;
  rayStack[Nrays] = r;
  ++Nrays;

  // keep looping until the stack is exhausted or the maximum number of rays is reached
  while(rayID<Nrays && Nrays<p_maxNrays){

    // get ray
    r = rayStack[rayID];
    
    // look for intersection of this ray with shapes
    int currentShapeID = -1;
    dfloat t = 20000.f;

    // look through grid to find intersections with ray
    gridRayIntersectionSearch(r, Nshapes, shapes, grid, &t, currentShapeID);
    
    // none found
    if(currentShapeID == -1){
      if(rayID==0)
	c = bg;

      // go to next ray
      ++rayID;
      continue;
    }

    // shape at nearest ray intersection
    shape_t currentShape = shapes[currentShapeID];
    
    // compute intersection location
    vector_t intersection = vectorAdd(r.start, vectorScale(t, r.dir));
    
    // find unit surface normal
    vector_t n = computeNormal(intersection, currentShape);
    
    /* use shadow tracing to determine color contribution from this intersection */
    dfloat rdotn = vectorDot(r.dir, n);

    /* Find the material to determine the colour */
    material_t currentMat = computeMaterial(Nmaterials, materials, intersection, currentShape);

    // test for reflection
    info_t info = currentMat.info;

    if(info.emitter==1){
      dfloat lambert = rdotn * r.coef; 
      c.red   += lambert * currentMat.diffuse.red;
      c.green += lambert * currentMat.diffuse.green;
      c.blue  += lambert * currentMat.diffuse.blue;
    }
    else{

      if(info.reflector==1){

	/* start ray slightly off surface */
	dfloat sc = p_shadowDelta;
	if(rdotn>0) // reverse offset if inside 
	  sc *= -1.f; // sign ? was -1
	
	vector_t shadowStart = vectorAdd(intersection, vectorScale(sc, n)); // HACK to shift ray start off service
	ray_t lightRay;
	lightRay.start = shadowStart;
	
	/* Find the value of the light at this point */
	for(unsigned int j=0; j < Nlights; j++){
	  
	  light_t currentLight = lights[j];
	  
	  vector_t dist = vectorSub(currentLight.pos, shadowStart);
	  if(vectorDot(n, dist) <= 0) continue;

	  dfloat lightDist = vectorNorm(dist);

	  dfloat tshadow = lightDist;
	  if(tshadow <= 0) continue;
	  
	  lightRay.dir = vectorScale((1.f/tshadow), dist);
	  
	  /* search in light ray direction for object */
	  int shadowShapeID = -1;
	  gridRayIntersectionSearch(lightRay, Nshapes, shapes, grid, &tshadow, shadowShapeID);

	  // check for objects in path of shadow ray
	  bool inShadow = false;	  
	  if(shadowShapeID==-1) // no object causes shadow
	    inShadow = false;
	  else if(tshadow >= 0 && tshadow < lightDist) // 
	    inShadow = true;

	  if(inShadow==false){
	    /* Lambert diffusion */
	    dfloat lambert = vectorDot(lightRay.dir, n) * r.coef;
	    c.red   += lambert * currentLight.intensity.red   * currentMat.diffuse.red;
	    c.green += lambert * currentLight.intensity.green * currentMat.diffuse.green;
	    c.blue  += lambert * currentLight.intensity.blue  * currentMat.diffuse.blue;
	  }
	}

	if((r.level+1<p_maxLevel) && Nrays<p_maxNrays) {
	  ray_t reflectRay;
	  // create new ray starting from offset intersection, with ray direction reflected in normal
	  reflectRay.start = shadowStart;
	  reflectRay.dir   = vectorAdd(r.dir, vectorScale(-2.0f*rdotn, n)); 

	  // increment level for new ray
	  reflectRay.level = r.level+1;
	  reflectRay.coef = r.coef*currentMat.reflection; // scale intensity
	
	  // launch new ray
	  rayStack[Nrays] = reflectRay;
	  // increment ray counter
	  ++Nrays;
	}
      }

      // https://www.scratchapixel.com/code.php?id=13&origin=/lessons/3d-basic-rendering/introduction-to-shading
      // test for refraction
      if(info.refractor==1){
	// can we add a new refraction ray to the stack ?
	if((r.level+1<p_maxLevel) && Nrays<p_maxNrays){
	
	  // push ray onto other side of surface
	  dfloat sc = -p_shadowDelta; // reverse number above
	  if(rdotn>0) 
	    sc *= -1;

	  // HACK to shift ray start off service
	  vector_t shadowStart = vectorAdd(intersection, vectorScale(sc, n)); 
	
	  // get index of refraction
	  dfloat eta = currentMat.eta; 
	
	  if(rdotn>0){
	    rdotn *= -1;
	  }else{
	    eta = 1.f/eta;
	  }
	
	  dfloat kappa = 1.f - eta*eta*(1.f - rdotn*rdotn);
	
	  if(kappa>0){
	    // create new refraction ray
	    ray_t refractRay;

	    // https://www.cs.cornell.edu/courses/cs4620/2012fa/lectures/36raytracing.pdf
	    // newdir = eta*d  - n*(eta*d.n - sqrt((1-eta*eta*(1-(d.n)^2))))
	    dfloat fac = eta*rdotn+sqrt(kappa);  // was - (NEED TO DOUBLE CHECK - other normal)
	    refractRay.start = shadowStart;
	    refractRay.dir = vectorNormalize(vectorSub(vectorScale(eta, r.dir), vectorScale(fac, n)));
	    refractRay.level = r.level+1;
	    refractRay.coef = r.coef*currentMat.refraction; // scale intensity
	    rayStack[Nrays] = refractRay;
	    ++Nrays;
	  }
	}
      }
    }
    
    // go to next ray on stack
    ++rayID;
  }

  return c;
  
}

__global__ void renderKernel (const int NI,
			      const int NJ,
			      const grid_t grid, 
			      const sensor_t sensor,
			      const int Nshapes,
			      const shape_t *shapes,
			      const int Nlights,
			      const light_t *lights,
			      const int Nmaterials,
			      const material_t *materials,
			      const dfloat costheta,
			      const dfloat sintheta,
			      const dfloat *randomNumbers,
			      unsigned char *img
			      ){
  
  const colour_t bg = sensor.bg;

  int I = threadIdx.x + blockDim.x*blockIdx.x;
  int J = threadIdx.y + blockDim.y*blockIdx.y;

  if(I<NI && J<NJ){
    ray_t r;
      
    dfloat coef = 1.0;
    int level = 0;

    // look at this: https://en.wikipedia.org/wiki/3D_projection

    // 2.5 location of sensor pixel
    colour_t c;

    //    dfloat randI = randomNumbers[I];
    //    dfloat randJ = randomNumbers[J];
    dfloat x0 = sensor.eyeX.x;
    dfloat y0 = sensor.eyeX.y;
    dfloat z0 = sensor.eyeX.z;

    // multiple rays emanating from sensor, passing through lens and focusing at the focal plane
    // 1. compute intersection of ray passing through lens center to focal plane

    // (sensorX + alpha*(lensC -sensorX)).sensorN = focalPlaneOffset
    // alpha = (focalOffset-s.sensorN)/( (lensC-s).sensorN) [ . dot product ]

    dfloat cx = BOXSIZE/2., cy =  BOXSIZE/2., cz = BOXSIZE/2;
    
    vector_t sensorN = vectorCrossProduct(sensor.Idir, sensor.Jdir);
    vector_t sensorX = sensorLocation(NI, NJ, I, J, sensor);
    dfloat   focalPlaneOffset = sensor.focalPlaneOffset;
    vector_t centralRayDir = vectorSub(sensor.lensC, sensorX);
    dfloat alpha = (focalPlaneOffset - vectorDot(sensorX, sensorN))/vectorDot(centralRayDir, sensorN);

    // 2. target
    vector_t targetX = vectorAdd(sensorX, vectorScale(alpha, centralRayDir));
    
    x0 = sensorX.x;
    y0 = sensorX.y;
    z0 = sensorX.z;
    
    // 3.  loop over vertical offsets on lens (thin lens)
    c.red = 0; c.green = 0; c.blue = 0;

    for(int samp=0;samp<p_Nsamples;++samp){

      // aperture width
      int sampId = (I+J*NI + samp*blockDim.x*blockDim.y)%NRANDOM;
      dfloat offI = randomNumbers[2*sampId+0]*p_apertureRadius;
      dfloat offJ = randomNumbers[2*sampId+1]*p_apertureRadius; 

      // choose random starting point on lens (assumes lens and sensor arre parallel)
      if(samp>0) { // primary ray
	x0 = sensor.lensC.x + offI*sensor.Idir.x + offJ*sensor.Jdir.x;
	y0 = sensor.lensC.y + offI*sensor.Idir.y + offJ*sensor.Jdir.y;
	z0 = sensor.lensC.z + offI*sensor.Idir.z + offJ*sensor.Jdir.z;
      }
      
      dfloat dx0 = targetX.x - x0;
      dfloat dy0 = targetX.y - y0;
      dfloat dz0 = targetX.z - z0;

      dfloat L0 = sqrt(dx0*dx0+dy0*dy0+dz0*dz0);
      dx0 = dx0/L0;
      dy0 = dy0/L0;
      dz0 = dz0/L0;
      
      r.start.x = costheta*(x0-cx) - sintheta*(z0-cz) + cx;
      r.start.y = y0;
      r.start.z = sintheta*(x0-cx) + costheta*(z0-cz) + cz;
      
      r.dir.x = costheta*dx0 - sintheta*dz0;
      r.dir.y = dy0;
      r.dir.z = sintheta*dx0 + costheta*dz0;

      colour_t newc =
	trace(grid, Nshapes, shapes, Nlights, lights, Nmaterials, materials, r, level, coef, bg);

      dfloat sc = (samp==0) ? p_primaryWeight: 1.f;
      c.red   += sc*newc.red;
      c.green += sc*newc.green;
      c.blue  += sc*newc.blue;
      
    }

    // primary weighted average
    c.red   /= (p_primaryWeight+p_Nsamples-1);
    c.green /= (p_primaryWeight+p_Nsamples-1);
    c.blue  /= (p_primaryWeight+p_Nsamples-1);

    // reverse vertical because of lensing
    img[(I + (NJ-1-J)*NI)*3 + 0] = (unsigned char)min(  c.red*255.0f, 255.0f);
    img[(I + (NJ-1-J)*NI)*3 + 1] = (unsigned char)min(c.green*255.0f, 255.0f);
    img[(I + (NJ-1-J)*NI)*3 + 2] = (unsigned char)min( c.blue*255.0f, 255.0f);
  }  
}


#define BLOCKSIZE 1024
#define LOGBLOCKSIZE 10

// https://en.wikipedia.org/wiki/Prefix_sum
// Hillis and Steele
// [ can be done with far fewer barriers ]
__global__ void startScanKernel(const int N,
				  const int *v,
				  int *scanv,
				  int *starts){

  __shared__ int s_v0[BLOCKSIZE];
  __shared__ int s_v1[BLOCKSIZE];

  int j = threadIdx.x;
  int b = blockIdx.x;
  int n = j + b*BLOCKSIZE;

  s_v0[j] = (n<N) ?  v[j+b*BLOCKSIZE]: 0;

  int offset = 1;
  do{
    __syncthreads();

    s_v1[j] = (j<offset) ? s_v0[j] : (s_v0[j]+s_v0[j-offset]) ;

    offset *= 2;

    __syncthreads();
    
    s_v0[j] = (j<offset) ? s_v1[j] : (s_v1[j]+s_v1[j-offset]) ;

    offset *= 2;
  } while(offset<BLOCKSIZE);

  if(n<N)
    scanv[n+1] = s_v0[j];

  if(j==(BLOCKSIZE-1)){
    starts[b+1] = s_v0[j];
  }
  
}

__global__ void finishScanKernel(const int N,
				 int *scanv,
				 int *starts){

  int j = threadIdx.x;
  int b = blockIdx.x;

  int n=j+b*BLOCKSIZE;

  if(n<N){
    int start = starts[b];
    
    scanv[n+1] += start;
  }
}

// returns the cumulative sum
int scan(const int N, const int *c_v, int *c_starts, int *starts, int *c_scanv){

  int B = BLOCKSIZE;
  int G = (N+BLOCKSIZE-1)/BLOCKSIZE;

  startScanKernel <<< G, B >>> (N, c_v, c_scanv, c_starts);

  cudaMemcpy(starts, c_starts, (G+1)*sizeof(int), cudaMemcpyDeviceToHost);

  starts[0] = 0;
  for(int b=0;b<G;++b){
    starts[b+1] += starts[b];
  }

  int count = starts[G];
  
  cudaMemcpy(c_starts, starts, (G+1)*sizeof(int), cudaMemcpyHostToDevice);

  finishScanKernel <<< G, B >>> (N, c_scanv, c_starts);

  return count;
}


__device__ bbox_t createBoundingBoxTriangle(triangle_t &triangle){

  bbox_t bbox;

  bbox.xmin = min(triangle.vertices[0].x, min(triangle.vertices[1].x, triangle.vertices[2].x));
  bbox.xmax = max(triangle.vertices[0].x, max(triangle.vertices[1].x, triangle.vertices[2].x));

  bbox.ymin = min(triangle.vertices[0].y, min(triangle.vertices[1].y, triangle.vertices[2].y));
  bbox.ymax = max(triangle.vertices[0].y, max(triangle.vertices[1].y, triangle.vertices[2].y));
  
  bbox.zmin = min(triangle.vertices[0].z, min(triangle.vertices[1].z, triangle.vertices[2].z));
  bbox.zmax = max(triangle.vertices[0].z, max(triangle.vertices[1].z, triangle.vertices[2].z));

  return bbox;
}

__device__ bbox_t createBoundingBoxSphere(sphere_t &sphere){

  bbox_t bbox;

  bbox.xmin = sphere.pos.x - sphere.radius;
  bbox.xmax = sphere.pos.x + sphere.radius;

  bbox.ymin = sphere.pos.y - sphere.radius;
  bbox.ymax = sphere.pos.y + sphere.radius;

  bbox.zmin = sphere.pos.z - sphere.radius;
  bbox.zmax = sphere.pos.z + sphere.radius;

  return bbox;
}

__device__ bbox_t createBoundingBoxEllipsoid(ellipsoid_t &ellipsoid){

  bbox_t bbox;

  bbox.xmin = ellipsoid.pos.x - 1.f/ellipsoid.invRadii.x;
  bbox.xmax = ellipsoid.pos.x + 1.f/ellipsoid.invRadii.x;

  bbox.ymin = ellipsoid.pos.y - 1.f/ellipsoid.invRadii.y;
  bbox.ymax = ellipsoid.pos.y + 1.f/ellipsoid.invRadii.y;

  bbox.zmin = ellipsoid.pos.z - 1.f/ellipsoid.invRadii.z;
  bbox.zmax = ellipsoid.pos.z + 1.f/ellipsoid.invRadii.z;

  return bbox;
}


__device__ bbox_t createBoundingBoxCylinder(cylinder_t &cylinder){

  bbox_t bbox;

  vector_t c = cylinder.center;
  vector_t a = cylinder.axis;
  dfloat   R = cylinder.radius;
  dfloat   H = cylinder.height;

  // xmax = c.x + (H/2)*a.x + (H/2)*|a.x| + R*|cross(e_x,a) |
  bbox.xmax = c.x + (H/2)*a.x + (H/2)*fabs(a.x) + R*sqrtf(a.y*a.y + a.z*a.z);
  bbox.xmin = c.x + (H/2)*a.x - (H/2)*fabs(a.x) - R*sqrtf(a.y*a.y + a.z*a.z);

  bbox.ymax = c.y + (H/2)*a.y + (H/2)*fabs(a.y) + R*sqrtf(a.x*a.x + a.z*a.z);
  bbox.ymin = c.y + (H/2)*a.y - (H/2)*fabs(a.y) - R*sqrtf(a.x*a.x + a.z*a.z);

  bbox.zmax = c.z + (H/2)*a.z + (H/2)*fabs(a.z) + R*sqrtf(a.x*a.x + a.y*a.y);
  bbox.zmin = c.z + (H/2)*a.z - (H/2)*fabs(a.z) - R*sqrtf(a.x*a.x + a.y*a.y);

  return bbox;
}

__device__ bbox_t createBoundingBoxCone(cone_t &cone){

  bbox_t bbox;

  vector_t v = cone.vertex;
  vector_t a = cone.axis;
  dfloat   R = cone.radius;
  dfloat   H = cone.height;

  bbox.xmax = max(v.x, v.x + H*a.x + R*sqrtf(a.y*a.y + a.z*a.z));
  bbox.xmin = min(v.x, v.x + H*a.x - R*sqrtf(a.y*a.y + a.z*a.z));

  bbox.ymax = max(v.y, v.y + H*a.y + R*sqrtf(a.x*a.x + a.z*a.z));
  bbox.ymin = min(v.y, v.y + H*a.y - R*sqrtf(a.x*a.x + a.z*a.z));

  bbox.zmax = max(v.z, v.z + H*a.z + R*sqrtf(a.x*a.x + a.y*a.y));
  bbox.zmin = min(v.z, v.z + H*a.z - R*sqrtf(a.x*a.x + a.y*a.y));
  
  return bbox;
}



__device__ bbox_t createBoundingBoxRectangle(rectangle_t &rectangle){

  bbox_t bbox;

  vector_t C = rectangle.center;
  vector_t A1 = rectangle.axis[0];
  vector_t A2 = rectangle.axis[1];
  vector_t n = vectorCrossProduct(A1, A2);

  dfloat   L1 = rectangle.length[0];
  dfloat   L2 = rectangle.length[1];
  A1 = vectorScale(L1/2., A1);
  A2 = vectorScale(L2/2., A2);
  
  dfloat delta = 1e-1;
  vector_t dn = vectorScale(delta, n);
  vector_t Cdown = vectorSub(C, dn);
  vector_t Cup = vectorAdd(C, dn);

  vector_t v[8];
  
  // C - delta*n + A1*(-L1/2) + A2*(-L2/2)
  v[0] = vectorSub(Cdown, vectorAdd(A1, A2));
  // C - delta*n + A1*(+L1/2) + A2*(-L2/2)
  v[1] = vectorAdd(Cdown, vectorSub(A1, A2));
  // C - delta*n + A1*(+L1/2) + A2*(+L2/2)
  v[2] = vectorAdd(Cdown, vectorAdd(A1, A2));
  // C - delta*n + A1*(-L1/2) + A2*(+L2/2)
  v[3] = vectorAdd(Cdown, vectorSub(A2, A1));

  // C + delta*n + A1*(-L1/2) + A2*(-L2/2)
  v[4] = vectorSub(Cup, vectorAdd(A1, A2));
  // C + delta*n + A1*(+L1/2) + A2*(-L2/2)
  v[5] = vectorAdd(Cup, vectorSub(A1, A2));
  // C + del6a*n + A1*(+L1/2) + A2*(+L2/2)
  v[6] = vectorAdd(Cup, vectorAdd(A1, A2));
  // C + delta*n + A1*(-L1/2) + A2*(+L2/2)
  v[7] = vectorAdd(Cup, vectorSub(A2, A1));

  bbox.xmin = 1e9;
  bbox.ymin = 1e9;
  bbox.zmin = 1e9;
  bbox.xmax = -1e9;
  bbox.ymax = -1e9;
  bbox.zmax = -1e9;

#pragma unroll 8
  for(int n=0;n<8;++n){
    bbox.xmin = min(bbox.xmin, v[n].x);
    bbox.ymin = min(bbox.ymin, v[n].y);
    bbox.zmin = min(bbox.zmin, v[n].z);
    bbox.xmax = max(bbox.xmax, v[n].x);
    bbox.ymax = max(bbox.ymax, v[n].y);
    bbox.zmax = max(bbox.zmax, v[n].z);
  }

  return bbox;
}

__device__ bbox_t createBoundingBoxDisk(disk_t &disk){

  bbox_t bbox;

  vector_t n = disk.normal;
  vector_t c = disk.center;
  dfloat   R = disk.radius;
  dfloat   H = .1; // assert thickness in normal

  // xmax = c.x + (H/2)*a.x + (H/2)*|a.x| + R*|cross(e_x,a) |
  bbox.xmax = c.x  + (H/2)*fabs(n.x) + R*sqrtf(n.y*n.y + n.z*n.z);
  bbox.xmin = c.x  - (H/2)*fabs(n.x) - R*sqrtf(n.y*n.y + n.z*n.z);

  bbox.ymax = c.y  + (H/2)*fabs(n.y) + R*sqrtf(n.x*n.x + n.z*n.z);
  bbox.ymin = c.y  - (H/2)*fabs(n.y) - R*sqrtf(n.x*n.x + n.z*n.z);

  bbox.zmax = c.z  + (H/2)*fabs(n.z) + R*sqrtf(n.x*n.x + n.y*n.y);
  bbox.zmin = c.z  - (H/2)*fabs(n.z) - R*sqrtf(n.x*n.x + n.y*n.y);

  return bbox;
}




__device__ void createBoundingBoxShape(const grid_t &grid, shape_t &shape){

  bbox_t bbox;

  switch(shape.type){
  case TRIANGLE:    bbox = createBoundingBoxTriangle(shape.triangle);   break;
  case SPHERE:      bbox = createBoundingBoxSphere(shape.sphere);       break;
  case ELLIPSOID:   bbox = createBoundingBoxEllipsoid(shape.ellipsoid); break;
  case IMAGE:
  case RECTANGLE:   bbox = createBoundingBoxRectangle(shape.rectangle); break;
  case CYLINDER:    bbox = createBoundingBoxCylinder(shape.cylinder);   break;
  case DISK:        bbox = createBoundingBoxDisk(shape.disk);           break;
  case CONE:        bbox = createBoundingBoxCone(shape.cone);           break;
  }

  int imin = floor(grid.invdx*(bbox.xmin-grid.xmin));
  int imax = floor(grid.invdx*(bbox.xmax-grid.xmin));
  int jmin = floor(grid.invdy*(bbox.ymin-grid.ymin));
  int jmax = floor(grid.invdy*(bbox.ymax-grid.ymin));
  int kmin = floor(grid.invdz*(bbox.zmin-grid.zmin));
  int kmax = floor(grid.invdz*(bbox.zmax-grid.zmin)); // was ceil
  
  bbox.imin = iclamp(imin, 0, grid.NI-1);
  bbox.imax = iclamp(imax, 0, grid.NI-1);

  bbox.jmin = iclamp(jmin, 0, grid.NJ-1);
  bbox.jmax = iclamp(jmax, 0, grid.NJ-1);

  bbox.kmin = iclamp(kmin, 0, grid.NK-1);
  bbox.kmax = iclamp(kmax, 0, grid.NK-1);
  
  shape.bbox = bbox;
  
}


__global__ void countShapesInBoxesKernel(const grid_t grid, const int Nshapes, shape_t *shapes, int *counts){

  int n = threadIdx.x + blockDim.x*blockIdx.x;

  if(n<Nshapes){
    shape_t &shape = shapes[n];
    createBoundingBoxShape(grid, shape);

    const  int imin = shape.bbox.imin;
    const  int imax = shape.bbox.imax;
    const  int jmin = shape.bbox.jmin;
    const  int jmax = shape.bbox.jmax;
    const  int kmin = shape.bbox.kmin;
    const  int kmax = shape.bbox.kmax;
    
    for(int k=kmin;k<=kmax;++k){
      for(int j=jmin;j<=jmax;++j){
	for(int i=imin;i<=imax;++i){
	  int id = i + j*grid.NI + k*grid.NI*grid.NJ;

	  atomicAdd(counts+id, 1);
	}
      }
    }
  }
}


__global__ void addShapesInBoxesKernel(const grid_t grid, const int Nshapes, const shape_t *shapes, int *boxCounters, int *boxContents){

  const int n = threadIdx.x + blockDim.x*blockIdx.x;
  
  if(n<Nshapes){

    const shape_t &shape = shapes[n];

    const  int imin = shape.bbox.imin;
    const  int imax = shape.bbox.imax;
    const  int jmin = shape.bbox.jmin;
    const  int jmax = shape.bbox.jmax;
    const  int kmin = shape.bbox.kmin;
    const  int kmax = shape.bbox.kmax;
    
    for(int k=kmin;k<=kmax;++k){
      for(int j=jmin;j<=jmax;++j){
	for(int i=imin;i<=imax;++i){
	  // box	  
	  const int id = i + j*grid.NI + k*grid.NI*grid.NJ;

	  // index in this box (post decremented)
	  // grab counter for this cell into index, then increment counter for this cell
	  const int index = atomicAdd(boxCounters+id,1); 

	  boxContents[index] = shape.id;
	}
      }
    }
  }
}


void populateGrid(grid_t *grid, int Nshapes, shape_t *c_shapes){

  if(grid->c_boxStarts){
    cudaFree(grid->c_boxStarts);
    cudaFree(grid->c_boxContents);
  }

  int Nboxes = grid->NI*grid->NJ*grid->NK;

  int *boxCounts = (int*) calloc(Nboxes+1, sizeof(int));
  
  int *c_boxCounts, *c_boxCounters;

  cudaMalloc(&c_boxCounts,         (Nboxes+1)*sizeof(int));
  cudaMalloc(&(grid->c_boxStarts), (Nboxes+1)*sizeof(int));
  cudaMalloc(&(c_boxCounters),     (Nboxes+1)*sizeof(int));    

  cudaMemset(c_boxCounts, 0, (Nboxes+1)*sizeof(int));

  int B = BLOCKSIZE;
  int G = (Nshapes+B-1)/B;

  countShapesInBoxesKernel <<< G, B >>> (*grid, Nshapes, c_shapes, c_boxCounts);
  
  // parallel scan to get cumulative offsets starting at zero
  int *tmp = (int*) calloc(Nboxes+1, sizeof(int));
  int *c_tmp;
  cudaMalloc(&c_tmp, (Nboxes+1)*sizeof(int));

  int Nentries =
    scan(Nboxes, c_boxCounts, c_tmp, tmp, grid->c_boxStarts);
  
  // build container for boxes
  cudaMalloc(&(grid->c_boxContents), (Nentries+1)*sizeof(int));

  cudaMemcpy(c_boxCounters, grid->c_boxStarts, (Nboxes+1)*sizeof(int), cudaMemcpyDeviceToDevice);
  
  // add each shape to every box that intersects the shape's bounding box
  addShapesInBoxesKernel <<< G, B >>> (*grid, Nshapes, c_shapes, c_boxCounters, grid->c_boxContents);
  
  free(tmp);
  free(boxCounts);

  cudaFree(c_tmp);
  cudaFree(c_boxCounts);
  cudaFree(c_boxCounters);
  
}

void sceneOffload(scene_t  *scene){

  grid_t *grid = scene->grid;
  
  cudaMalloc(&(scene->c_materials), scene->Nmaterials*sizeof(material_t));
  cudaMalloc(&(scene->c_lights),    scene->Nlights*sizeof(light_t));
  cudaMalloc(&(scene->c_shapes),    scene->Nshapes*sizeof(shape_t));

  cudaMalloc(&(scene->c_img), WIDTH*HEIGHT*3*sizeof(char));

  cudaMalloc(&(grid->c_bboxes), (grid->NI*grid->NJ*grid->NK)*sizeof(bbox_t));

  cudaMemcpy(scene->c_shapes,    scene->shapes,    scene->Nshapes*sizeof(shape_t),       cudaMemcpyHostToDevice);
  cudaMemcpy(scene->c_materials, scene->materials, scene->Nmaterials*sizeof(material_t), cudaMemcpyHostToDevice);
  cudaMemcpy(scene->c_lights,    scene->lights,    scene->Nlights*sizeof(light_t),       cudaMemcpyHostToDevice);

  cudaMemcpy(grid->c_bboxes,      grid->bboxes, (grid->NI*grid->NJ*grid->NK)*sizeof(bbox_t), cudaMemcpyHostToDevice);

  cudaMalloc(&(scene->c_randomNumbers), 2*NRANDOM*sizeof(dfloat));

  cudaMemcpy(scene->c_randomNumbers, scene->randomNumbers, 2*NRANDOM*sizeof(dfloat), cudaMemcpyHostToDevice);

  
}


dfloat drandRange48(dfloat dmin, dfloat dmax){

  return dmin + drand48()*(dmax-dmin);
}


// L = size of box
// delta = width of layer around box
grid_t *gridSetup(dfloat L, dfloat delta){ 
  
  // bu	ild grid
  grid_t *grid = (grid_t*) calloc(1, sizeof(grid_t));

  grid->xmin = -delta;
  grid->xmax = L + delta;
  grid->ymin = -delta;
  grid->ymax = L + delta;
  grid->zmin = -delta;
  grid->zmax = L + delta;
  
  grid->NI = 401;
  grid->NJ = 401;
  grid->NK = 401;

  grid->dx = (grid->xmax-grid->xmin)/grid->NI;
  grid->dy = (grid->ymax-grid->ymin)/grid->NJ;
  grid->dz = (grid->zmax-grid->zmin)/grid->NK;
  
  grid->invdx = grid->NI/(grid->xmax-grid->xmin);
  grid->invdy = grid->NJ/(grid->ymax-grid->ymin);
  grid->invdz = grid->NK/(grid->zmax-grid->zmin);

  grid->bboxes = (bbox_t*) calloc(grid->NI*grid->NJ*grid->NK, sizeof(bbox_t));
  for(int k=0;k<grid->NK;++k){
    for(int j=0;j<grid->NJ;++j){
      for(int i=0;i<grid->NI;++i){
	int id = i + j*grid->NI + k*grid->NI*grid->NJ;
	grid->bboxes[id].xmin = i*grid->dx + grid->xmin;
	grid->bboxes[id].xmax = (i+1)*grid->dx + grid->xmin;
	grid->bboxes[id].ymin = j*grid->dy + grid->ymin;
	grid->bboxes[id].ymax = (j+1)*grid->dy + grid->ymin;
	grid->bboxes[id].zmin = k*grid->dz + grid->zmin;
	grid->bboxes[id].zmax = (k+1)*grid->dz + grid->zmin;
      }
    }
  }

  return grid;
}



scene_t *sceneSetup(int  plotNelements,
		    dfloat *plotx,
		    dfloat *ploty,
		    dfloat *plotz,
		    dfloat *plotq){

  int i;
  
  int Nmaterials = 64;
  material_t *materials = (material_t*) calloc(Nmaterials, sizeof(material_t));

  materials[0].diffuse.red   = 1;
  materials[0].diffuse.green = 1;
  materials[0].diffuse.blue  = 1;
  materials[0].reflection = 1;
  materials[0].eta = 1;
  materials[0].refraction = 0;

  materials[0].info.refractor = 0;
  materials[0].info.reflector = 1;
  materials[0].info.emitter = 0;

  materials[1].diffuse.red = 0;
  materials[1].diffuse.green = 240/255.;
  materials[1].diffuse.blue = 20/255.;
  materials[1].reflection = .3;
  materials[1].eta = .7;
  materials[1].refraction = .1;
  materials[1].info.refractor = 1;
  materials[1].info.reflector = 1;
  materials[1].info.emitter = 0;


  for(i=2;i<Nmaterials;++i){
    dfloat red = 0,green = 0,blue =0;

    red   = drandRange48(0.125,0.8);
    green = drandRange48(0.125,0.8);
    blue  = drandRange48(0.125,0.8);
    
    materials[i].diffuse.red   = red;
    materials[i].diffuse.green = green;
    materials[i].diffuse.blue  = blue;

    materials[i].eta = 2;
    materials[i].refraction = 1; // transmission coeff
    if(drand48() > .5){
      materials[i].reflection = .9;
      materials[i].info.reflector = 1;
    }
    if(drand48() > .5){
      materials[i].refraction = .9;
      materials[i].info.refractor = 1;
    }

    if(!materials[i].info.refractor && !materials[i].info.reflector){
      materials[i].info.reflector = 1;
    }

#if 0
    printf("materials[%d] = {{rgb=%g,%g,%g},{reflection=%g,refraction=%g,eta=%g},info {reflector=%d,refractor=%d,emitter=%d}\n",
	   i,
	   materials[i].diffuse.red,
	   materials[i].diffuse.green,
	   materials[i].diffuse.blue,
	   materials[i].reflection,
	   materials[i].refraction,
	   materials[i].eta,
	   materials[i].info.reflector,
	   materials[i].info.refractor,
	   materials[i].info.emitter);
#endif    
  }

  int Ntriangles = plotNelements;
  
  int Nrectangles = 1;
  int Nshapes = Ntriangles + Nrectangles;

  // length of side of world box
  dfloat L = BOXSIZE;
  
  shape_t *shapes = (shape_t*) calloc(Nshapes, sizeof(shape_t));

  dfloat triXmin = 1e9, triXmax = -1e9;
  dfloat triYmin = 1e9, triYmax = -1e9;	
  dfloat triZmin = 1e9, triZmax = -1e9;

  for(int n=0;n<plotNelements*3;++n){
    triXmin = min(triXmin, plotx[n]);
    triXmax = max(triXmax, plotx[n]);
    triYmin = min(triYmin, ploty[n]);
    triYmax = max(triYmax, ploty[n]);
    triZmin = min(triZmin, plotz[n]);
    triZmax = max(triZmax, plotz[n]);
  }

  printf("Ntriangles = %d in range (%lg,%lg x %lg,%lg x %lg,%lg)\n",
	 Ntriangles, triXmin, triXmax, triYmin, triYmax, triZmin, triZmax);
  
  dfloat maxL = max(triXmax-triXmin, max(triYmax-triYmin, triZmax-triZmin));
  
  int bcnt = 0;

  dfloat brot  = 0;
  dfloat bcosrot = cos(brot);
  dfloat bsinrot = sin(brot);
  dfloat boffx = 500; drandRange48(250, L-250);
  dfloat boffy = 500; drandRange48(250, L-250);
  dfloat boffz = 500; drandRange48(250, L-250);
  dfloat bscal = 800;
  int bmat = 32; 

  printf("bmat = %d\n", bmat);

  dfloat newTriXmin = 1e9, newTriXmax = -1e9;
  dfloat newTriYmin = 1e9, newTriYmax = -1e9;	
  dfloat newTriZmin = 1e9, newTriZmax = -1e9;

  for(i=0;i<Ntriangles;++i){

    for(int v=0;v<3;++v){
      shapes[bcnt].triangle.vertices[v] =
	vectorCreate((plotx[i*3+v]-triXmin)/maxL,
		     (plotz[i*3+v]-triZmin)/maxL,
		     (ploty[i*3+v]-triYmin)/maxL); // swapped y and z

      shapes[bcnt].triangle.q[v] = plotq[i*3+v];
    }
    vector_t tmp = shapes[bcnt].triangle.vertices[1];
    shapes[bcnt].triangle.vertices[1] = shapes[bcnt].triangle.vertices[2];
    shapes[bcnt].triangle.vertices[2] = tmp;
    
    for(int v=0;v<3;++v){
      
      dfloat x = bscal*shapes[bcnt].triangle.vertices[v].x;
      dfloat y = L - bscal*shapes[bcnt].triangle.vertices[v].y;
      dfloat z = bscal*shapes[bcnt].triangle.vertices[v].z;
      
      dfloat xrot =  cos(brot)*x + sin(brot)*z;
      dfloat zrot = -sin(brot)*x + cos(brot)*z;
      
      shapes[bcnt].triangle.vertices[v].x = boffx + xrot;
      shapes[bcnt].triangle.vertices[v].y = y;
      shapes[bcnt].triangle.vertices[v].z = boffz + zrot;

#if 1
      newTriXmin = min(newTriXmin, shapes[bcnt].triangle.vertices[v].x);
      newTriXmax = max(newTriXmax, shapes[bcnt].triangle.vertices[v].x);
      newTriYmin = min(newTriYmin, shapes[bcnt].triangle.vertices[v].y);
      newTriYmax = max(newTriYmax, shapes[bcnt].triangle.vertices[v].y);
      newTriZmin = min(newTriZmin, shapes[bcnt].triangle.vertices[v].z);
      newTriZmax = max(newTriZmax, shapes[bcnt].triangle.vertices[v].z);
#endif
    }
    
    shapes[bcnt].material = bmat;
    shapes[bcnt].type = TRIANGLE;
    shapes[bcnt].id = bcnt;
    ++bcnt;
  }


  printf("Ntriangles = %d in range (%lg,%lg x %lg,%lg x %lg,%lg)\n",
	 Ntriangles, newTriXmin, newTriXmax, newTriYmin, newTriYmax, newTriZmin, newTriZmax);
  
  int cnt = Ntriangles;
  
  // add one rectangle
  if(Nrectangles>0){

    vector_t a = vectorCreate(0, L, 0);
    vector_t b = vectorCreate(0, L, L);
    vector_t c = vectorCreate(L, L, L);
    vector_t d = vectorCreate(L, L, 0);

    vector_t ab = vectorSub(d,a);
    vector_t ad = vectorSub(b,a);
    shapes[cnt].rectangle.length[0] = vectorNorm(ab);
    shapes[cnt].rectangle.length[1] = vectorNorm(ad);
    shapes[cnt].rectangle.axis[0] = vectorNormalize(ab);
    shapes[cnt].rectangle.axis[1] = vectorNormalize(ad);
    shapes[cnt].rectangle.center  = vectorScale(0.25, vectorAdd(vectorAdd(a,b),vectorAdd(c,d)));
    shapes[cnt].material = -1;
    shapes[cnt].type = RECTANGLE;
    shapes[cnt].id = cnt;    
    ++cnt;
  }

  int Nlights = 5;
  light_t *lights = (light_t*) calloc(Nlights, sizeof(light_t));
	
  lights[0].pos.x = L/2;
  lights[0].pos.y = 0;
  lights[0].pos.z = -100;
  lights[0].intensity.red = 1;
  lights[0].intensity.green = 1;
  lights[0].intensity.blue = 1;
	
  lights[1].pos.x = 3200;
  lights[1].pos.y = 3000;
  lights[1].pos.z = -1000;
  lights[1].intensity.red = 0.6;
  lights[1].intensity.green = 0.7;
  lights[1].intensity.blue = 1;

  lights[2].pos.x = 600;
  lights[2].pos.y = 0;
  lights[2].pos.z = -100;
  lights[2].intensity.red = 0.3;
  lights[2].intensity.green = 0.5;
  lights[2].intensity.blue = 1;

  lights[3].pos.x = L/2;
  lights[3].pos.y = 0;
  lights[3].pos.z = L/2;
  lights[3].intensity.red = 0.8;
  lights[3].intensity.green = 0.8;
  lights[3].intensity.blue = 1;

  lights[4].pos.x = L;
  lights[4].pos.y = L;
  lights[4].pos.z = -1000;
  lights[4].intensity.red = 1;
  lights[4].intensity.green = 1;
  lights[4].intensity.blue = 1;

  scene_t *scene = (scene_t*) calloc(1, sizeof(scene_t));

  scene->Ntriangles = plotNelements;
  
  scene->Nlights  = Nlights;
  scene->lights   = lights;
  scene->Nshapes   = Nshapes;
  scene->shapes   = shapes;
  scene->Nmaterials = Nmaterials;
  scene->materials  = materials;
  scene->grid = gridSetup(L, 1600);

  scene->randomNumbers = (dfloat*) calloc(2*NRANDOM, sizeof(dfloat));
  for(int i=0;i<NRANDOM;++i){
    dfloat r1 = 2*drand48()-1;
    dfloat r2 = 2*drand48()-1;

    scene->randomNumbers[2*i+0] = r1/sqrt(r1*r1+r2*r2);
    scene->randomNumbers[2*i+1] = r2/sqrt(r1*r1+r2*r2);
  }

  
  return scene;
}

/* Output data as PPM file */
void saveppm(char *filename, unsigned char *img, int width, int height){

  /* FILE pointer */
  FILE *f;
  
  /* Open file for writing */
  f = fopen(filename, "wb");
  
  /* PPM header info, including the size of the image */
  fprintf(f, "P6 %d %d %d\n", width, height, 255);

  /* Write the image data to the file - remember 3 byte per pixel */
  fwrite(img, 3, width*height, f);

  /* Make sure you close the file */
  fclose(f);
}


scene_t *simpleRayTracerSetup(int     plotNelements,
			      dfloat *plotx,
			      dfloat *ploty,
			      dfloat *plotz,
			      dfloat *plotq){

  
  // initialize triangles and spheres
  scene_t *scene = sceneSetup(plotNelements, plotx, ploty, plotz, plotq);
  
  // port to GPU
  sceneOffload(scene);
  
  int TX = 8, TY = 8;
  dim3 B(TX,TY,1);
  dim3 G( (WIDTH+TX-1)/TX, (HEIGHT+TY-1)/TY, 1);
  
  void populateGrid(grid_t *grid, int Nshapes, shape_t *c_shapes);
  populateGrid(scene->grid, scene->Nshapes, scene->c_shapes);

  return scene;
}




// to compile animation:
//   ffmpeg -y -i image_%05d.ppm -pix_fmt yuv420p foo.mp4

scene_t *scene = NULL;

void simpleRayTracer(int     plotNelements,
		     dfloat *plotx,
		     dfloat *ploty,
		     dfloat *plotz,
		     dfloat *plotq,
		     const char *fileBaseName,
		     const int fileIndex){
  
  // initialize triangles and spheres
  if(!scene)
    scene = simpleRayTracerSetup(plotNelements, plotx, ploty, plotz, plotq);
  
  // update field
  
  for(int i=0;i<scene->Ntriangles;++i){
    for(int v=0;v<3;++v){
      scene->shapes[i].triangle.q[v] = plotq[i*3+v];
    }
  }
  
  cudaMemcpy(scene->c_shapes,    scene->shapes,    scene->Nshapes*sizeof(shape_t),       cudaMemcpyHostToDevice);
  
  // 1. location of observer eye (before rotation)
  sensor_t sensor;
  
  // background color
  sensor.bg.red   = 126./256;
  sensor.bg.green = 192./256;
  sensor.bg.blue  = 238./256;
  
  dfloat br = 3.75f*BOXSIZE;
  
  // angle elevation to y-z plane
  dfloat eyeAngle = .5*M_PI/2.f; // 0 is above, pi/2 is from side.  M_PI/3; 0; M_PI/2.;
  
  // target view
  vector_t targetX = vectorCreate(BOXSIZE/2., BOXSIZE, BOXSIZE/2.); // this I do not understand why target -B/2
  sensor.eyeX = vectorAdd(targetX, vectorCreate(0, -br*cos(eyeAngle), -br*sin(eyeAngle))); 
  dfloat sensorAngle = eyeAngle; +15.*M_PI/180.;
  sensor.Idir   = vectorCreate(1.f, 0.f, 0.f);
  sensor.Jdir   = vectorCreate(0.f, sin(sensorAngle), -cos(sensorAngle));
  vector_t sensorNormal = vectorCrossProduct(sensor.Idir, sensor.Jdir);
  
#if 0
  printf("eyeX = %g,%g,%g \n, IDir = %g,%g,%g, \n, Jdir = %g,%g,%g, Ndir = %g,%g,%g\n",
	 sensor.eyeX.x,
	 sensor.eyeX.y,
	 sensor.eyeX.z,
	 sensor.Idir.x,
	 sensor.Idir.y,
	 sensor.Idir.z,
	 sensor.Jdir.x,
	 sensor.Jdir.y,
	 sensor.Jdir.z,
	 sensorNormal.x,
	 sensorNormal.y,
	 sensorNormal.z);
#endif	 
    
  // 2.4 length of sensor in axis 1 & 2
  sensor.Ilength = 20.f;
  sensor.Jlength = HEIGHT*20.f/WIDTH;
  sensor.offset  = 0.f;
  
  // 2.5 normal distance from sensor to focal plane
  dfloat lensOffset = 50;
  sensor.lensC = vectorAdd(sensor.eyeX, vectorScale(lensOffset, vectorCrossProduct(sensor.Idir, sensor.Jdir)));
  
  // why 0.25 ?
  sensor.focalPlaneOffset = 0.22f*fabs(vectorTripleProduct(sensor.Idir, sensor.Jdir, vectorSub(targetX,sensor.eyeX))); // triple product
    
  //  printf("lensOffset = %g, sensor.focalPlaneOffset = %g\n", lensOffset, sensor.focalPlaneOffset);
    
  /* rotation angle in y-z */
  dfloat theta = M_PI*fileIndex*1./180.;
    
  int TX = 8, TY = 8;
  dim3 B(TX,TY,1);
  dim3 G( (WIDTH+TX-1)/TX, (HEIGHT+TY-1)/TY, 1);
    
  /* render scene */
  renderKernel
    <<<G,B>>>(WIDTH,
	      HEIGHT,
	      scene->grid[0],
	      sensor,
	      scene->Nshapes,
	      scene->c_shapes,
	      scene->Nlights,
	      scene->c_lights,
	      scene->Nmaterials,
	      scene->c_materials,
	      cos(theta), 
	      sin(theta),
	      scene->c_randomNumbers,
	      scene->c_img);
    
  /* copy image back to host */
  unsigned char *img = (unsigned char*) calloc(3*WIDTH*HEIGHT, sizeof(char));
  cudaMemcpy(img, scene->c_img, 3*WIDTH*HEIGHT*sizeof(char), cudaMemcpyDeviceToHost);
  
  // make sure images directory exists
  mkdir("images", S_IRUSR | S_IREAD | S_IWUSR | S_IWRITE | S_IXUSR | S_IEXEC);
  
  char fileName[BUFSIZ];
  
  sprintf(fileName, "images/%s_%05d.ppm", fileBaseName, fileIndex);
  saveppm(fileName, img, WIDTH, HEIGHT);
  
  free(img);
  
}
