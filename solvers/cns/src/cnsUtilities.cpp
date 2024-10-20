#include "cns.hpp"
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

int countinString(std::string s, char c){
    // Count variable
    int res = 0;
    for (int i=0;i< (int)s.length();i++){
        if (s[i] == c){res++;}
    }
    return res;
}



void cns_t::setFlowStates(){
  char comma        = ',' ;  // seperates the entries
  char semicolon    = ';' ;  // seperates the groups
  // Read Referenece State
  std::string str, sword, dword;
  settings.getSetting("FLOW STATES", str);

  std::stringstream ss(str);
  std::stringstream ssc(str);
  // Count comma seperated reference state data
  NstatePoints = (2+mesh.dim);   // Density - Velocity - Pressure
  NstateSets   = countinString(ss.str(), semicolon);

  while(ssc){
     getline(ssc, sword, semicolon);
     if((int) sword.length()>0){
     int Ndata = countinString(sword, comma)+1;
     LIBP_ABORT("correct the number of inputs in tokenizer", NstatePoints!=Ndata);
    }
  }
  ssc.clear(); 


  flowStates.malloc(NstatePoints*NstateSets, 0.0); 
  // Check the correct number of data points
  int sks = 0, skd = 0; 
  while(!ss.eof()){
    getline(ss, sword, semicolon);
    if((int)sword.length()>0){
    std::stringstream ds(sword);
     while(!ds.eof()){
      getline(ds, dword, comma);
      if(dword.length()>0){
        flowStates[sks*NstatePoints + skd] = std::stod(dword); 
        // printf("%.4e \n", flowStates[sks*NstatePoints + skd]); 
        skd++; 
      }
     }
     skd=0; sks++; 
   }
  }

  ss.clear();

  o_flowStates = platform.malloc<dfloat>(flowStates);
  // IC State State ID
  settings.getSetting("IC STATE ID", ICStateID); 
  // BC State ID
  settings.getSetting("BC STATE ID", BCStateID); 
  // REF State ID
  settings.getSetting("REFERENCE STATE ID", RefStateID); 
}


void cns_t::setBoundaryMaps(){
  char comma        = ',' ;  // seperates the entries
  char colon        = ':' ;  // seperates the entries
  char semicolon    = ';' ;  // seperates the entries
  // Read Referenece State
  std::string str, sword, tword, dword;
  settings.getSetting("GEOMETRIC-TO-PHYSICAL MAP", str);

  std::stringstream ss(str);
  std::stringstream ssc(str);
  
  // Number of geoemtric sets
  // int NgeoSets   = countinString(ss.str(), colon);
  NgeoIDs = 0; 
  while(!ssc.eof()){
     getline(ssc, sword, semicolon);
     if((int)sword.length()>0){
        NgeoIDs += countinString(sword, comma)+1;
      }
  }

  GToB.malloc(NgeoIDs*2, 0);

   int bcid=0; int skg = 0; 
  // Check the correct number of data points
  while(!ss.eof()){
    getline(ss, sword, semicolon);
    if((int)sword.length()>0){    
    int sk1=0;
    std::stringstream ts(sword); 
     while(!ts.eof()){
      getline(ts, tword, colon);
      // std::cout<< "tword: "<<tword <<" " << std::endl; 
      if((int)tword.length()>0){
        if(sk1==0){
          bcid = std::stoi(tword); 
        }else{
         std::stringstream ds(tword);
          while(!ds.eof()){
            getline(ds, dword, comma);
             if((int)tword.length()>0){
               GToB[skg*2 + 0] = std::stoi(dword);
               GToB[skg*2 + 1] = bcid;
               // std::cout<< GToB[skg*2 + 0]<< " "<<GToB[skg*2 + 1]<<std::endl; 
               skg++;
             }
           }
        }
        sk1++;
      }
     }
   }
  }
  ss.clear(); 
  ssc.clear();


  // Convert BCType to 
  EToG.malloc(mesh.Nelements*mesh.Nfaces,0); 
  for(int e=0; e<mesh.Nelements; e++){
    for(int f=0; f<mesh.Nfaces; f++){
      const int id = e*mesh.Nfaces + f; 
      int gtype = mesh.EToB[id]; 
      EToG[id]  = gtype; 
      if(gtype>0){
        for(int i=0; i<NgeoIDs; i++){
          int gmapID = GToB[2*i +0]; 
          int bmapID = GToB[2*i +1];
          if(gtype==gmapID){
            mesh.EToB[id] = bmapID; 
          }
        }
      }
    }
  }
  mesh.o_EToB.copyFrom(mesh.EToB); 

  o_EToG = platform.malloc<int>(EToG);

  // Number of reference points per state
  props["defines/" "p_NsP"]    = (int) NstatePoints;
  props["defines/" "p_NsS"]    = (int) NstateSets;

}



void cns_t::setReport(){
  char comma        = ',' ;  // seperates the entries
  char colon        = ':' ;  // seperates the entries
  char semicolon    = ';' ;  // seperates the entries

  reportForces     = settings.compareSetting("REPORT FORCES", "TRUE") ? 1:0;
  reportMoments    = settings.compareSetting("REPORT MOMENTS","TRUE") ? 1:0;
  reportComponent  = settings.compareSetting("REPORT COMPONENT","TRUE") ? 1:0;

  // Read Referenece State
  std::string str, sword, tword, dword;
  settings.getSetting("REPORT GROUPS", str);

  std::stringstream ss(str);
  std::stringstream ssc(str);
  
  // Number of geoemtric sets
  NreportGroups   = countinString(ss.str(), colon);  
  NreportIDs = 0; 
  while(!ssc.eof()){
     getline(ssc, sword, semicolon);
     if(sword.length()>0){
        NreportIDs += countinString(sword, comma)+1;
      }
  }
  reportGroups.malloc(2*NreportIDs, 0);

   int gid=0; int skg = 0; 
  // Check the correct number of data points
  while(!ss.eof()){
    getline(ss, sword, semicolon);
    if((int)sword.length()>0){    
    int sk1=0;
    std::stringstream ts(sword); 
     while(!ts.eof()){
      getline(ts, tword, colon);
      if((int)tword.length()>0){
        if(sk1==0){
          gid = std::stoi(tword); 
        }else{
         std::stringstream ds(tword);
          while(!ds.eof()){
            getline(ds, dword, comma);
             if((int)tword.length()>0){
               reportGroups[skg*2 + 0] = std::stoi(dword);
               reportGroups[skg*2 + 1] = gid;
               // std::cout<< NreportGroups<< " "<<NreportIDs<< " "<< reportGroups[skg*2 + 0]<< " "<<reportGroups[skg*2 + 1]<<std::endl; 
               skg++;
             }
           }
        }
        sk1++;
      }
     }
   }
  }

  o_reportGroups = platform.malloc<int>(reportGroups);



  // Read Referenece State
  settings.getSetting("MOMENT CENTER", str);
  std::stringstream ssm(str);
  std::stringstream ssmc(str); 
  
  // printf("%s \n", str.c_str());
  // Number of geoemtric sets
  int Nmoment   = countinString(ssm.str(), comma);  
  LIBP_ABORT("correct the number of inputs in tokenizer", Nmoment!=mesh.dim);
  momentCenter.malloc(2*Nmoment, 0);
  
  int sk=0; 
  while(ssmc){
    getline(ssmc, sword, comma);
     if((int)sword.length()>0){   
      momentCenter[sk++] =  std::stod(sword);
      // printf("%s %.4e \n", sword.c_str(), momentCenter[sk-1]);
    }
  }
  // LIBP_ABORT("correct the number of inputs in tokenizer", sk!=mesh.dim);

  o_momentCenter = platform.malloc<dfloat>(momentCenter);

  ss.clear(); ssc.clear();
  ssm.clear(); ssmc.clear();
  // Number of reference points per state
  props["defines/" "p_NrGrp"]    = (int) NreportGroups;
  props["defines/" "p_NrIDs"]    = (int) NreportIDs;

}