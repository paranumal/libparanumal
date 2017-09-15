#include <iostream>

#include "occa.hpp"

int main(int argc, char **argv){
  //occa::printAvailableDevices();

  int N = atoi(argv[1]);
  int Kblk = atoi(argv[2]);
  //  printf("N = %d, Kblk = %d\n",N,Kblk);
  int Np = (N+1)*(N+2)*(N+3)/6;

  int Nq;
  switch (N){
  case 1: Nq = 4; //Nq = 6;
    break;
  case 2: Nq = 11; //Nq = 14;
    break;
  case 3: Nq = 23; //Nq = 31;
    break;
  case 4: Nq = 44; //Nq = 57;
    break;
  case 5: Nq = 74; //Nq = 95;
    break;
  case 6: Nq = 122; //Nq = 146;
    break;
  case 7: Nq = 177; //Nq = 214;
    break;
  }
  int K = 5e4;

  float *Vq = new float[Nq*Np];
  float *Pq = new float[Nq*Np];
  for (int i = 0; i < Np*Nq; ++i){
    Vq[i] = i/Np;
    Pq[i] = i/Nq;
  }

  float *wq = new float[Nq*K];
  for (int i = 0; i < Nq*K; ++i){
    wq[i] = 1.f;
  }

  float *PK = new float[Np*Np*K];
  for (int i = 0; i < Np*Np*K; ++i){
    PK[i] = 1.f;
  }

  float *u = new float[Np*K];
  for (int i = 0; i < Np*K; ++i){
    u[i] = (float) i;
  }

  occa::device device;
  occa::kernelInfo info;
  info.addDefine("p_Nfields", 1);
  info.addDefine("p_Np", Np);
  info.addDefine("p_Np2", Np*Np);
  info.addDefine("p_cubNp", Nq);
  info.addDefine("p_Kblk", Kblk);

  //---[ Device setup with string flags ]-------------------
  //device.setup("mode = Serial");
  // device.setup("mode = OpenMP  , schedule = compact, chunk = 10");
  // device.setup("mode = OpenCL  , platformID = 0, deviceID = 1");
  device.setup("mode = CUDA    , deviceID = 2");

  occa::memory o_u, o_wq,  o_Vq, o_Pq, o_PK;
  o_u  = device.malloc(Np*K*sizeof(float),u);
  o_wq = device.malloc(Nq*K*sizeof(float),wq);
  o_Vq = device.malloc(Nq*Np*sizeof(float),Vq);
  o_Pq = device.malloc(Nq*Np*sizeof(float),Pq);
  o_PK = device.malloc(Np*Np*K*sizeof(float),PK);

  // OKL: OCCA Kernel Language
  occa::kernel applyWADG = device.buildKernelFromSource("wadgKernels.okl","applyWADG",info);
  occa::kernel applyPK = device.buildKernelFromSource("wadgKernels.okl","applyPw",info);

  int nsteps = 100;
  using namespace std;
  clock_t begin = clock();
  for (int i = 0; i < nsteps; ++i){
    applyWADG(K, o_Vq, o_Pq, o_wq, o_u);
  }
  device.finish();
  clock_t end = clock();
  double wadg_time = double(end - begin) / (CLOCKS_PER_SEC);

  begin = clock();
  for (int i = 0; i < nsteps; ++i){
    applyPK(K, o_PK, o_u);
  }
  device.finish();
  end = clock();
  double PK_time = double(end - begin) / (CLOCKS_PER_SEC);

  double scale = nsteps*K; // compute time per element
  printf("N = %d, Kblk = %d: wadg time = %g seconds, PK time = %g seconds\n",N,Kblk,wadg_time/scale,PK_time/scale);

  device.free();

  return 0;
}
