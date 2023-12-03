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

#include "cns.hpp"


void myTokenizer(const int NrefState, std::string s, memory<dfloat> & refState, char del)
{
    std::stringstream ss(s);
    std::string word;
    int nref = 0; 
    while (!ss.eof()) {
        getline(ss, word, del);
        refState[nref] = std::stod(word); 
        nref ++; 
    }

    LIBP_ABORT("Correct the reference state: pref, uref, vref, wref, pref", nref!= NrefState);
  
}


void cns_t::setupPhysics(properties_t & props){
 
  // Set isentropic exponent and related info
  settings.getSetting("GAMMA", gamma);
  
  // These are needed to simplfy computations 
  gammaM1  = gamma - 1.0; 
  gammaP1  = gamma + 1.0; 
  igamma   = 1.0/ gamma; 
  igammaM1 = 1.0/ gammaM1; 
  igammaP1 = 1.0/ gammaP1; 

  // Read Referenece State
  std::string stateStr;
  const int NrefState = 6;
  refState.malloc(NrefState, 0.0); 
  settings.getSetting("REFERENCE STATE", stateStr);
  myTokenizer(NrefState, stateStr,  refState,  ',');

  // dfloat rRef =0.0, vRef   = 0.0, pRef = 0.0, tRef = 0.0, mRef = 0.0, lRef = 0.0; 
  // rRef = refState[0]; 
  // vRef = refState[1]; 
  // pRef = refState[2]; 
  // tRef = refState[3]; 
  // mRef = refState[4]; 
  // lRef = refState[5]; 


  // Set specific gas constant
  if(settings.compareSetting("NONDIMENSIONAL EQUATIONS", "TRUE")){
      settings.getSetting("REYNOLDS NUMBER", Re);    
      settings.getSetting("MACH NUMBER", Ma); 
      R = 1.0 / (gamma*Ma*Ma); 
  }else {
      settings.getSetting("SPECIFIC GAS CONSTANT", R); 
  }

  // Set presure and volumetric expansion coefficients
  cp = R*gamma/(gamma-1.0);  
  cv = R/(gamma-1.0); 

  if(settings.compareSetting("SOLVER TYPE", "NAVIER-STOKES")){
    settings.getSetting("PRANDTL NUMBER", Pr); 
    settings.getSetting("VISCOSITY", mu);
    // update mu as 1/Re
    mu = (settings.compareSetting("NONDIMENSIONAL EQUATIONS", "TRUE"))? 1.0/Re : mu;
  }else{
    mu = 0.0; Re = 0.0; 
  }

  Nph  = 8;  
  pCoeff.malloc(Nph,0.0);
  MUID = 0; 
  GMID = 1; 
  PRID = 2; 
  RRID = 3; 
  CPID = 4; 
  CVID = 5; 
  KAID = 6; 
  M2ID = 7; 
  pCoeff[MUID] = mu; // Bulk Viscosity
  pCoeff[PRID] = Pr; // Prandtl Number
  pCoeff[RRID] = R;  // Specific Gas Constant
  pCoeff[GMID] = gamma; 
  pCoeff[CPID] = cp;
  pCoeff[CVID] = cv;
  pCoeff[KAID] = 1.0/( (gamma-1.0)*Ma*Ma)*mu/Pr;
  pCoeff[M2ID] = Ma*Ma;


  if(settings.compareSetting("SOLVER TYPE", "NAVIER-STOKES")){
    if(settings.compareSetting("VISCOSITY TYPE", "CONSTANT")){
      viscType = 1;
    }else if(settings.compareSetting("VISCOSITY TYPE", "SUTHERLAND")){
      viscType  = 2;
      Nph      += 4;  
      EXID = 6; 
      TRID = 7; 
      TSID = 8; 
      CSID = 9; 
      pCoeff.realloc(Nph);
      // Coefficients from White, F. M., Viscous fluid flow, McGraw-Hill, 2006
      dfloat Tref = 273.15, Ts = 110.4, exp = 1.5;  
      pCoeff[EXID] = exp;            // exponent  
      pCoeff[TRID] = 1.0/Tref;       // inverse of reference temperature here !!!    
      pCoeff[TSID] = Ts/Tref;        // Ts/Tref approximately !!! 
      pCoeff[CSID] = pow(pCoeff[TSID],pCoeff[EXID])*(1.0+pCoeff[TSID])/(2.0*pCoeff[TSID]*pCoeff[TSID]); // exponent  
    }else if(settings.compareSetting("VISCOSITY TYPE", "POWER-LAW")){
      viscType  = 3;
      dfloat exp = 1.5;
      dfloat Tref = (settings.compareSetting("NONDIMENSIONAL EQUATIONS", "TRUE")) ? 
                    1.0 : 273.15; 
      pCoeff[MUID] = mu / Tref * exp; // Update viscosity
      pCoeff[EXID] = exp;             // exponent  
      pCoeff[TRID] = Tref;            // inverse of reference temperature = 1/Tref   
    }
  }else{ // Euler solver
    viscType = 0; 
  }

  

  // Define physical model on Device 
  props["defines/" "p_viscType"]= viscType;
  props["defines/" "p_MUID"]    = MUID;
  props["defines/" "p_GMID"]    = GMID;
  props["defines/" "p_PRID"]    = PRID;
  props["defines/" "p_RRID"]    = RRID;
  props["defines/" "p_CPID"]    = CPID;
  props["defines/" "p_CVID"]    = CVID;
  props["defines/" "p_KAID"]    = KAID;
  props["defines/" "p_M2ID"]    = M2ID;
  props["defines/" "p_EXID"]    = EXID;
  props["defines/" "p_TRID"]    = TRID;
  props["defines/" "p_TSID"]    = TSID;
  props["defines/" "p_CSID"]    = CSID;

  // Define reference state on device
  // props["defines/" "p_RBAR"]    = refState[0];
  // props["defines/" "p_UBAR"]    = refState[1];
  // props["defines/" "p_VBAR"]    = refState[2];
  // if(mesh.dim==3){
  //   props["defines/" "p_WBAR"]    = refState[3];
  //   props["defines/" "p_PBAR"]    = refState[4];

  // }else{
  //   props["defines/" "p_PBAR"]    = refState[3];    
  // }
 
  // move physical model to device
  o_pCoeff = platform.malloc<dfloat>(pCoeff);   

}
