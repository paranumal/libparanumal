/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#ifndef CNS_HPP
#define CNS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"
// #include "stab.hpp"

#define DCNS LIBP_DIR"/solvers/cns/"

using namespace libp;

class cnsSettings_t: public settings_t {
public:
  cnsSettings_t(comm_t& _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);

  // stabSettings_t extractStabSettings();
};



namespace Stab {
  enum Detect { NODETECT =0, ALL =1, KLOCKNER =2, PERSSON  =3, DUCROS   =4,};
  enum Type { NOSTAB  =0, FILTER  =1, LIMITER =2, ARTDIFF =3, SUBCELL =4,};
} //namespace stab



class cns_t: public solver_t {
public:
  mesh_t mesh;
  // stab_t stab {};
  // stabSettings_t stabSettings; 

  properties_t props; 

  int Nfields;
  int Ngrads;

 

  dfloat Re, Ma, Pr;
  dfloat mu, cp, cv, R, kappa, gamma;
  // dfloat gamma, igamma, gammaM1, igammaM1, gammaP1, igammaP1;

  dfloat beta_ldg, tau_ldg;

  int viscType;
  int Nph; // number of physical parameters 
  int MUID, GMID, RRID, PRID, CPID; 
  int CVID, KAID, M2ID, BTID, TAID; 
  int EXID, TRID, TSID, CSID;  

   // Number of reference points
  int NstatePoints, NstateSets, NgeoIDs;  
  int ICStateID, BCStateID;  
  // Physical coefficients
  memory<dfloat> flowStates; 
  deviceMemory<dfloat> o_flowStates;

  // Physical coefficients
  memory<dfloat> GToB; 
  memory<int>   EToG;  
  deviceMemory<int> o_EToG; 
    
  // moment center
  int reportForces, reportMoments; 
  int reportComponent, NreportGroups, NreportIDs; 
  memory<dfloat> momentCenter; 
  deviceMemory<dfloat> o_momentCenter;
  memory<int>  reportGroups;      
  deviceMemory<int>  o_reportGroups;  

  


  // Physical coefficients
  memory<dfloat> pCoeff; 
  deviceMemory<dfloat> o_pCoeff; 

  int cubature;
  int isothermal;

  timeStepper_t timeStepper;

  ogs::halo_t fieldTraceHalo;
  ogs::halo_t gradTraceHalo;

  memory<dfloat> q;
  deviceMemory<dfloat> o_q;

  kernel_t volumeKernel;
  kernel_t surfaceKernel;
  kernel_t cubatureVolumeKernel;
  kernel_t cubatureSurfaceKernel;

  kernel_t gradVolumeKernel;
  kernel_t gradSurfaceKernel;

  kernel_t cubatureGradVolumeKernel;
  kernel_t cubatureGradSurfaceKernel;

  kernel_t vorticityKernel;

  kernel_t constrainKernel;

  kernel_t initialConditionKernel;
  kernel_t maxWaveSpeedKernel;

  kernel_t limiterReconstructKernel; 
  kernel_t limiterGradientKernel; 
  kernel_t limiterVertexBoundaryKernel;

  kernel_t forcesVolumeKernel; 
  kernel_t forcesSurfaceKernel; 

  kernel_t reportArrangeLayoutKernel; 
  kernel_t reportAverageKernel; 

  // for post processing only 
  ogs::ogs_t ogs;
  memory<dfloat> weight; 
  deviceMemory<dfloat> o_weight; 

  // Stabilization Related
  int detectType=0, stabType=0; 

  // Detected elements and detected field variable
  memory<dlong>elmList; 
  deviceMemory<dlong> o_elmList; 
  deviceMemory<int> o_modeIds; 

  // always scaled such that 0<=s<=1.0 
  memory<dfloat> qdetect; 
  deviceMemory<dfloat>o_qdetect; 
  deviceMemory<dfloat>o_invV, o_LS1D, o_BLD; 

  memory<dfloat>  projectC0, qducros; 
  deviceMemory<dfloat>o_projectC0, o_qducros; 


  // Artificial Viscosity
  dfloat avisAlpha; 
  memory<dfloat> elmLength, Vmax; 
  deviceMemory<dfloat> o_elmLength, o_Vmax;


  memory<dfloat> vertexViscosity, viscosity; 
  deviceMemory<dfloat> o_vertexViscosity, o_viscosity;

  memory<dfloat>PM12N;
  deviceMemory<dfloat>o_PM12N;


  kernel_t detectKernel;
  kernel_t detectDucrosKernel;
  kernel_t computeViscosityKernel;
  kernel_t projectViscosityKernel;
  kernel_t maxVelocityKernel;


  // Artificial viscosity 


 
  cns_t() = default;
  cns_t(platform_t &_platform, mesh_t &_mesh,
              cnsSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             cnsSettings_t& _settings);

  void Run();
  // Set reference values, thermodynamics, nondimensional values
  void setupPhysics();

  // void tokenizer(const int N, std::string s, memory<dfloat> & state, char delimiter);
  // template<typename T>
  
  void writeForces(dfloat time, int tstep, int frame);
  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  
  // Setup detectors and stabilizers
  void setupStab();
  void applyStab(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_gradQ, const dfloat T);
  void applyDetect(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_gradQ, const dfloat T);

  void setupKlocknerDetector();
  void applyKlocknerDetector(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_gradQ, const dfloat T);

  void rhsNoStab(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
  void rhsArtDiff(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
  void rhsLimiter(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);



  void setupArtificialDiffusion();
  void applyArtificialDiffusion(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_gradQ, const dfloat T);
  void setupLimiter();
  void setupNoStab();
  
  void Report(dfloat time, int tstep) override;

  void PlotFields(memory<dfloat> Q, memory<dfloat> V, std::string fileName);

  dfloat MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T);

  // Utilities
  void setFlowStates();
  void setBoundaryMaps();
  void setReport();

dfloat ElementViscosityScaleTri2D(dlong e); 
dfloat ElementViscosityScaleQuad2D(dlong e); 
dfloat ElementViscosityScaleTet3D(dlong e); 
dfloat ElementViscosityScaleHex3D(dlong e); 

};

#endif
