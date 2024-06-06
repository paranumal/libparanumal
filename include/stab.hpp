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

#ifndef STAB_HPP
#define STAB_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "settings.hpp"
#include "mesh.hpp"
#include "ogs.hpp"

namespace libp{

class stabSettings_t: public settings_t {
public:
  stabSettings_t() = default;
  stabSettings_t(const comm_t& _comm);
  void report();
};


namespace Stab {
  /*Element types*/
  enum Solver {
    HJS    = 1,
    CNS    = 2,
    INS    = 3,
  };

  enum Detector {
    NODETECT =0,
    ALL      =1,
    KLOCKNER =2,
    PERSSON  =3,
    DUCROS   =4,
  };

  enum Type {
    NOSTAB  =0, 
    FILTER  =1,
    LIMITER =2,
    ARTDIFF =3,
    SUBCELL =4,
  };

} //namespace Stab

class stab_t {
 public:
  platform_t platform;
  properties_t props;
  stabSettings_t settings;
  mesh_t mesh; 
  comm_t comm;

  ogs::halo_t traceHalo;

  // # of detector and stabiization fields;
  int Ndfields, Nsfields;  

  /***************************/
  int CXID, CYID, CZID, IVID; 
  int FXID, FYID, FZID, SAID; 
  int NXID, NYID, NZID, BCID; 
  /***************************/
  int VID; 

  /*************************/
  /* Solver Data           */
  /*************************/
  Stab::Solver solver;
  Stab::Detector detector;
  Stab::Type type;

  /*************************/
  /* Detector Data        */
  /*************************/
  memory<dlong> elementList;
  deviceMemory<dlong> o_elementList; 

  memory<dfloat> elementListF;
  deviceMemory<dfloat> o_elementListF; 

  memory<dfloat> qdetector;
  deviceMemory<dfloat> o_qdetector; 

  memory<dfloat> qducros;
  deviceMemory<dfloat> o_qducros; 
  
  //     Klockner Detector      // 
  memory<int> modeMap; 
  deviceMemory<int> o_modeMap;
  memory<dfloat> leastSquares1D; 
  memory<dfloat> baseLineDecay; 

  deviceMemory<dfloat> o_leastSquares1D; 
  deviceMemory<dfloat> o_baseLineDecay;
  deviceMemory<dfloat> o_invV;
 
  // Persson Detector
  memory<dfloat> projectNm1; 
  deviceMemory<dfloat> o_projectNm1; 

  kernel_t copyToIntKernel; 
  kernel_t copyToFloatKernel; 
  kernel_t detectKernel; 
  kernel_t extractFieldKernel; 

  kernel_t detectKernel2; 

  /*****************************/
  /*   FILTER STABILIZATION    */
  /*****************************/
  memory<dfloat> filterM; 
  deviceMemory<dfloat> o_filterM; 

  // kernel_t filterKernel; 
 


  memory<dfloat> weight; 
  deviceMemory<dfloat> o_weight; 
  /*****************************/
  /*   LIMITER STABILIZATION    */
  /*****************************/
  memory<dfloat> qv; 
  deviceMemory<dfloat> o_qv; 

  memory<dfloat> qc; 
  deviceMemory<dfloat> o_qc; 

  memory<dfloat> dq; 
  deviceMemory<dfloat> o_dq; 

  memory<dfloat> dqf; 
  deviceMemory<dfloat> o_dqf; 

  // Project to cell centers
  memory<dfloat> projectC0; 
  deviceMemory<dfloat> o_projectC0; 

  memory<int> vertexNodes;
  deviceMemory<int> o_vertexNodes; 

  memory<dlong> lvmapM, lvmapP; 
  deviceMemory<dfloat> o_lvmapM, o_lvmapP; 


  kernel_t projectCellAverageKernel;
  kernel_t getVertexValuesKernel;

  /*******************************************/
  /*   ARTIFICIAL DIFFUSION STABILIZATION    */
  /*******************************************/
  dfloat scaleFactor;  
  // Memory for artificial Viscosity Activation Function 
  memory<dfloat> viscosityScale; 
  deviceMemory<dfloat> o_viscosityScale; 

  

  dlong Nmasked; 
  memory<dlong> elementMask; 
  deviceMemory<dlong> o_elementMask; 

  // Memory for artificial Viscosity Activation Function 
  memory<dfloat> viscosityActivation; 
  deviceMemory<dfloat> o_viscosityActivation; 

  memory<dfloat> vertexViscosity; 
  deviceMemory<dfloat> o_vertexViscosity; 

  memory<dfloat> viscosity; 
  deviceMemory<dfloat> o_viscosity; 

  memory<dfloat> projectViscosity; 
  deviceMemory<dfloat> o_projectViscosity; 

  // smmoth out viscosity
  ogs::ogs_t ogs;
  
  // Vertex Viscosity  
  kernel_t computeViscosityKernel; 
  kernel_t projectViscosityKernel; 
  kernel_t maskElementsKernel; 

  /*******************************************/
  /*         SUBCELL STABILIZATION           */
  /*******************************************/
  int N, Nsubcells, Nint, Next, Nfields; 
  int Nverts, Nfaces, NfaceVertices, Np; 
  int Nvgeo, Nsgeo; 

  // minor girid connectivity
  memory<int> mEToV, mEToE, mEToF, mFaceNodes; 
  memory<int> faceVertices; 
  memory<dfloat> vr, vs, vt; // local vertex coordinates @ reference element
  memory<dfloat> cr, cs, ct; // center points of subcells @ reference element 
  memory<dfloat> fr, fs, ft; // center points of subcells @ reference element 
  memory<dfloat> mJ; // Jacobian of minor grid

  memory<dfloat> sq, sqf, srhsq; 
  deviceMemory<dfloat> o_sq, o_sqf, o_srhsq;

  memory<int> mFToE, mFToF, mDGID; 
  deviceMemory<int> o_mFToE, o_mFToF, o_mDGID, o_mEToV, o_mFaceNodes;   

  // local projection 
  memory<dfloat> PM,  RM, PVM; // volume reconstuction, projection,  vertex  projection
  deviceMemory<dfloat> o_PM, o_RM, o_PVM;  

  memory<dfloat> PFM, RFM, SLIFT; // face projection and reconstruction 
  deviceMemory<dfloat> o_PFM, o_RFM, o_SLIFT; 

  // Connectivity
  memory<int> ielist, eelist; // local connectivity 
  deviceMemory<int> o_ielist, o_eelist; 

  // Glocal Connectivity
  memory<dlong> emapP, fmapP; // global connectivity
  deviceMemory<dlong> o_emapP, o_fmapP; // global connectivity

  // geometric info for limiter and subcell
  memory<dfloat> vgeo, sgeo; 
  deviceMemory<dfloat> o_vgeo, o_sgeo; 

  deviceMemory<dlong> o_EToE; 

  kernel_t findNeighKernel; 
  kernel_t projectFVKernel;
  kernel_t reconstructDGKernel;


  stab_t() = default;
  stab_t(platform_t& _platform, mesh_t &_mesh, stabSettings_t& _settings, properties_t& _kernelInfo){
    Setup(_platform, _mesh, _settings, _kernelInfo);
  }

  // stab setup
  void Setup(platform_t& _platform, mesh_t &_mesh, stabSettings_t& _settings, properties_t& kernelInfo);

  void detectSetup(){
    switch (detector) {
      case Stab::KLOCKNER:
         detectSetupKlockner(); 
        break;
      case Stab::PERSSON:
         detectSetupPersson(); 
        break;
      case Stab::DUCROS:
         if(solver!=Stab::CNS){
          LIBP_FORCE_ABORT("DUCROS indicator is defined for CNS only");
         }
         detectSetupDucros(); 
        break;
      case Stab::ALL:
         detectSetupAll(); 
        break;
      case Stab::NODETECT:
         detectSetupNodetect(); 
        break;
    } 
  }

  void Detect(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
    switch (detector) {
      case Stab::KLOCKNER:
         detectApplyKlockner(o_Q, o_RHS, T); 
        break;
      case Stab::PERSSON:
         detectApplyPersson(o_Q, o_RHS, T); 
        break;
      case Stab::DUCROS:
         detectApplyDucros(o_Q, o_RHS, T); 
        break;
      case Stab::ALL:
         detectApplyAll(o_Q, o_RHS, T); 
        break;
      case Stab::NODETECT:
         detectApplyNodetect(o_Q, o_RHS, T); 
        break;
    } 
  }



  // It is safer o split solver implementations I guess
  void stabSetup(){
    switch (type) {
      case Stab::NOSTAB:
        break; 
      case Stab::FILTER:
        // stabSetupFilter();
        break;
      case Stab::LIMITER:
         stabSetupLimiter();
        break;
      case Stab::ARTDIFF:
        stabSetupArtdiff();
        break;
      case Stab::SUBCELL:
        // stabSetupSubcell();
        break;
    } 
  }



  // It is safer o split solver implementations I guess
  void Apply(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
    switch (type) {
      case Stab::NOSTAB:
        break;  
      case Stab::FILTER:
        // stabApplyFilter();
        break;
      case Stab::LIMITER:
         stabApplyLimiter(o_Q, o_RHS, T);
         // LIBP_FORCE_ABORT("Limiter is not implemented yet");
        break;
      case Stab::ARTDIFF:
        stabApplyArtdiff(o_Q, o_RHS, T);
        break;
      case Stab::SUBCELL:
        // stabApplySubcell();
        break;
    } 
  }



   //***********************GENERAL IO OPERATIONS****************************//
   void Report(dfloat time, int tstep); 
   void Test();



   dlong GetElementNumber(deviceMemory<dlong>& eList); 
   void PlotElements(memory<dlong> eList, const std::string fileName);
   void PlotFields(memory<dfloat> Q, const std::string fileName);

   // *******************Detector Related Functions***************************//
   void detectSetupAll(); 
   void detectApplyAll(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void detectSetupNodetect(); 
   void detectApplyNodetect(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void detectSetupKlockner(); 
   void detectApplyKlockner(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void detectSetupPersson(); 
   void detectApplyPersson(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void detectSetupDucros(); 
   void detectApplyDucros(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void ModeInfoKlocknerTri2D(int _N, memory<int>& _modeMap);
   void ModeInfoKlocknerQuad2D(int _N, memory<int>& _modeMap);
   void ModeInfoKlocknerTet3D(int _N, memory<int>& _modeMap);
   void ModeInfoKlocknerHex3D(int _N, memory<int>& _modeMap);

   void ModeInfoPerssonTri2D(int _N, memory<dfloat>&  _truncModes);
   void ModeInfoPerssonQuad2D(int _N, memory<dfloat>& _truncModes);
   void ModeInfoPerssonTet3D(int _N, memory<dfloat>&  _truncModes);
   void ModeInfoPerssonHex3D(int _N, memory<dfloat>&  _truncModes);

   void LeastSquaresFitKlockner(int _N, memory<dfloat>& _LSF);
   void BaseLineDecayKlockner(int _N, memory<dfloat>& _BLD);

   // *********************FILTER RELATED*****************************************//
   void stabSetupFilter(); 
   void stabApplyFilter(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void FilterMatrixTri2D (int _N, int _Nc, int _s, memory<dfloat>& _filterMatrix);
   void FilterMatrix1D(int _N, int _Nc, int _s, memory<dfloat>& _filterMatrix);
   void FilterMatrixTet3D (int _N, int _Nc, int _s, memory<dfloat>& _filterMatrix);

   // *********************ARTIFICIAL DIFFUSION RELATED*****************************//
   void stabSetupArtdiff(); 
   void stabApplyArtdiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 
    
   dfloat ElementViscosityScaleTri2D(dlong e);
   dfloat ElementViscosityScaleQuad2D(dlong e);
   dfloat ElementViscosityScaleTet3D(dlong e);
   dfloat ElementViscosityScaleHex3D(dlong e);
// *********************LIMITER RELATED*****************************//
   void stabSetupLimiter(); 
   void limiterSetupOperators(); 
   void stabApplyLimiter(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 




   void limiterGeometricFactorsTri2D(); 
   void limiterGeometricFactorsQuad2D(); 
   void limiterGeometricFactorsTet3D(); 
   void limiterGeometricFactorsHex3D(); 

   // *********************SUBCELL RELATED*****************************//
   void stabSetupSubcell(); 
   void stabSetupSubcellTri2D(); 
   void stabApplySubcell(); 

          // ******************************************************************//
   void CellLocalConnectTri2D(); 
   void CellGlobalConnectTri2D(); 
   void CellCreateMinorGridTri2D();
   void CellGeometricFactorsTri2D(); 
   void CellSetupOperatorsTri2D(); 
   void CellFindBestMatchTri2D(dfloat x1, dfloat y1, dlong eP, int fP, 
                               memory<int> &elist, memory<dfloat> &x2,  memory<dfloat> &y2, 
                               int &nE, int &nF);


   void CellEquispacedEToVTri2D(const int _N, memory<int>& mEToV);
   void CellWarpBlendEToVTri2D(const int _N, memory<int>& mEToV);

  
 private:
  /*Set the type of mesh*/
  void setTypes(const Stab::Solver solver, 
                const Stab::Detector detector, 
                const Stab::Type type);

 
};



}
#endif

