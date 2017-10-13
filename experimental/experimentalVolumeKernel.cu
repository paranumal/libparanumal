#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

// CUDA experimental cubature nonlinear term

#define dfloat float
#define p_N 5
#define p_Np ((int)((p_N+1)*(p_N+2))/2) 
#define p_cubNp ((int)(3*p_Np))
#define p_NSIMD 32
#define p_BSIMD ((int)((p_Np+p_NSIMD-1)/p_NSIMD))
#define p_CSIMD ((int)((p_cubNp+p_NSIMD-1)/p_NSIMD))

#define p_Nvgeo 7
#define p_RXID 0
#define p_RYID 1
#define p_SXID 2
#define p_SYID 3

__global__ void experimentalVolumeKernel(const int Nelements, 
					 const dfloat * __restrict__ vgeo,
					 const dfloat * __restrict__ cI,
					 const dfloat * __restrict__ cDr,
					 const dfloat * __restrict__ cDs,
					 const dfloat * __restrict__ cProj,
					 const dfloat * __restrict__ u,
					 const dfloat * __restrict__ v,
					 dfloat * __restrict__ Nu,
					 dfloat * __restrict__ Nv){
  
  const unsigned int mask  = 0xFFFFFFFF;

  const int e = blockIdx.x; 
  const int t = threadIdx.x;
  dfloat r_u[p_BSIMD];
  dfloat r_v[p_BSIMD];
  dfloat r_Nu[p_BSIMD];
  dfloat r_Nv[p_BSIMD];
  
  for(int b=0;b<p_BSIMD;++b){ // loop over blocks of 32
    const int n = t + b*p_NSIMD;
    const int id = e*p_Np+n;
    // load q once and stash in register
    r_u[b] = (n<p_Np) ? u[id] : 0.0f;
    r_v[b] = (n<p_Np) ? v[id] : 0.0f;
    r_Nu[b] = 0;
    r_Nv[b] = 0;
  }

  for(int c=0;c<p_CSIMD;++c){ // for each cubature node on this thread
    
    // compute u,v,dudr,duds,dvdr,dvds at each cubature node in this block then use them immediately
    
    dfloat cubu = 0;
    dfloat cubv = 0;
    dfloat cubur = 0;
    dfloat cubus = 0;
    dfloat cubvr = 0;
    dfloat cubvs = 0;
    
    const int i = t + c*p_NSIMD;

#if 1
    for(int b=0;b<p_BSIMD;++b){ // loop over blocks of 32
      
      for(int s=0;s<p_NSIMD;++s){
	
	const int m = s + b*p_NSIMD;
	const dfloat um = __shfl_sync(mask, r_u[b], s, p_NSIMD);
	const dfloat vm = __shfl_sync(mask, r_v[b], s, p_NSIMD);
	
	if(i<p_cubNp){
	  const dfloat cIim  =  cI[i+m*p_Np]; // weak -- needs L1
	  const dfloat cDrim = cDr[i+m*p_Np];
	  const dfloat cDsim = cDs[i+m*p_Np];
	  cubu  +=  cIim*um; 
	  cubur += cDrim*um; 
	  cubus += cDsim*um; 
	  
	  cubv  +=  cIim*vm; 
	  cubvr += cDrim*vm; 
	  cubvs += cDsim*vm; 
	}
      }
    }
#endif
    
    const dfloat rx = vgeo[e*p_Nvgeo + p_RXID];
    const dfloat ry = vgeo[e*p_Nvgeo + p_RYID];
    const dfloat sx = vgeo[e*p_Nvgeo + p_SXID];
    const dfloat sy = vgeo[e*p_Nvgeo + p_SYID];
    
    // now have u,v,ur,us,vr,vs at cubature node c
    const dfloat cubux = rx*cubur + sx*cubus;
    const dfloat cubuy = ry*cubur + sy*cubus;
    const dfloat cubvx = rx*cubvr + sx*cubvs;
    const dfloat cubvy = ry*cubvr + sy*cubvs;
    const dfloat cubNu = cubu*cubux + cubv*cubuy;
    const dfloat cubNv = cubu*cubvx + cubv*cubvy;

#if 1
    for(int b=0;b<p_BSIMD;++b){ 
      
      for(int s=0;s<p_NSIMD;++s){
	
	// not sure about this part
	const int m = s + c*p_NSIMD; 
	const int i = t + b*p_NSIMD;
	       
	const dfloat Num = __shfl_sync(mask, cubNu, s, p_NSIMD);
	const dfloat Nvm = __shfl_sync(mask, cubNv, s, p_NSIMD);

	if(i<p_Np && m<p_cubNp){
	  const dfloat cPRim = cProj[i+m*p_Np]; // weak - needs L1
	  r_Nu[b]  +=  cPRim*Num; 
	  r_Nv[b]  +=  cPRim*Nvm;
	}
      }
    }
#endif
  }
 
#if 0
  for(int b=0;b<p_BSIMD;++b){ // loop over blocks of 32
    const int id = t + b*p_NSIMD;
    if(id<p_Np){
      Nu[id + e*p_Np] = r_Nu[b];
      Nv[id + e*p_Np] = r_Nv[b];
    }
  }
#endif

}

unsigned long long Nbytes = 0;

void randAlloc(int Nrand, dfloat **h_v, dfloat **c_v){
  
  *h_v = (dfloat*) calloc(Nrand, sizeof(dfloat));
  for(int n=0;n<Nrand;++n) h_v[0][n] = drand48();
  
  cudaMalloc(c_v, Nrand*sizeof(dfloat));
  
  cudaMemcpy(*c_v, *h_v, Nrand*sizeof(dfloat), cudaMemcpyHostToDevice);

  Nbytes += Nrand*sizeof(dfloat);
  
}


int main(int argc, char **argv){

  int Nelements = atoi(argv[1]);

  dfloat *h_u, *h_v, *h_Nu, *h_Nv, *h_vgeo, *h_cI, *h_cDr, *h_cDs, *h_cProj ;
  dfloat *c_u, *c_v, *c_Nu, *c_Nv, *c_vgeo, *c_cI, *c_cDr, *c_cDs, *c_cProj ;

  randAlloc(p_Np*Nelements, &h_u, &c_u);
  randAlloc(p_Np*Nelements, &h_v, &c_v);

  randAlloc(p_Np*Nelements, &h_Nu, &c_Nu);
  randAlloc(p_Np*Nelements, &h_Nv, &c_Nv);

  randAlloc(p_Np*p_cubNp, &h_cI, &c_cI);

  randAlloc(p_Np*p_cubNp, &h_cDr, &c_cDr);
  randAlloc(p_Np*p_cubNp, &h_cDs, &c_cDs);

  randAlloc(p_Np*p_cubNp, &h_cProj, &c_cProj);

  randAlloc(p_Nvgeo*Nelements, &h_vgeo, &c_vgeo);

  printf("Nbytes = %llu\n", Nbytes);
  
  dim3 G(Nelements,1,1);
  dim3 B(p_NSIMD,1,1);

  printf("p_Np = %d\n", p_Np);
  printf("p_cubNp = %d\n", p_cubNp);
  printf("p_NSIMD = %d\n", p_NSIMD);
  printf("p_CSIMD = %d\n", p_CSIMD);
  printf("p_BSIMD = %d\n", p_BSIMD);
  printf("G.x = %d, B.x = %d\n", G.x, B.x);
  
  
  experimentalVolumeKernel <<< G, B >>> (Nelements, c_vgeo, c_cI, c_cDr, c_cDs, c_cProj, c_u, c_v, c_Nu, c_Nv);

  exit(0);
  return 0;
  
}
  
