
/*
  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void printEntry(FILE *fp, const char *name, const char *val){
  fprintf(fp, "[%s]\n%s\n\n", name, val);
}

void printEntry(FILE *fp, const char *name, const int val){
  fprintf(fp, "[%s]\n%d\n\n", name, val);
}

void printEntry(FILE *fp, const char *name, const double val){
  fprintf(fp, "[%s]\n%lf\n\n", name, val);
}


void printSetup(const char *device, int degree, int nx, int ny, int nz, double epsy){

  FILE *fp = fopen("setupHex3D.rc", "w");

  char *thread_model=strdup(device);
  
  char rcformat[]="2.0";
  char data_file[]="../solvers/elliptic/data/ellipticSine3D.h";
  char mesh[]="BOX";
  int dim=3, element=4, boundary_flag=1;
  char map_file[]="../solvers/elliptic/data/kershaw.okl";
  double map_param_y=epsy;
  int  map_model=1;

  int platform_number=0, device_number=0;
  double Lambda=0.0;
  char discretization[]="CONTINUOUS";
  char discretization_quadrature[]="GLL";
  char linear_solver[]="PCG";
  char precon[]="SEMFEM";
  char multigrid_smoother[]="CHEBYSHEV";
  int multigrid_cheby_degree=1;
  char multigrid_coarsening[]="HALFDEGREES";
  char paralmond_device_matrix_type[]="float";  
  char paralmond_cycle[]="VCYCLE";
  char paralmond_strength[]="RUGESTUBEN";
  double paralmond_strength_cutoff=0.5;
  char paralmond_aggregation[]="SMOOTHED";
  char paralmond_smoother[]="CHEBYSHEV";
  int paralmond_cheby_degree=2;
  char output_to_file[]="FALSE";
  
  printEntry(fp,"FORMAT", "2.0");
  printEntry(fp,"DATA FILE", data_file);
  printEntry(fp,"MESH FILE", "BOX");
  printEntry(fp,"MESH DIMENSION", 3);
  printEntry(fp,"ELEMENT TYPE",   12);
  printEntry(fp,"BOX NX", nx);
  printEntry(fp,"BOX NY", ny);
  printEntry(fp,"BOX NZ", nz);
  printEntry(fp,"BOX DIMX", 1);
  printEntry(fp,"BOX DIMY", 1);
  printEntry(fp,"BOX DIMZ", 1);
  printEntry(fp,"BOX BOUNDARY FLAG", 1);
  printEntry(fp,"BOX COORDINATE MAP FILE", map_file);
  printEntry(fp,"BOX COORDINATE MAP PARAMETER Y", map_param_y);
  printEntry(fp,"BOX COORDINATE MAP PARAMETER Z", map_param_y);
  printEntry(fp,"BOX COORDINATE MAP MODEL", map_model);
  printEntry(fp,"POLYNOMIAL DEGREE", degree);
  printEntry(fp,"THREAD MODEL", thread_model);
  printEntry(fp,"PLATFORM NUMBER", platform_number);
  printEntry(fp,"DEVICE NUMBER", device_number);
  printEntry(fp,"DISCRETIZATION", discretization);
  printEntry(fp,"ELLIPTIC INTEGRATION", discretization_quadrature);
  printEntry(fp,"LINEAR SOLVER", linear_solver);
  printEntry(fp,"LAMBDA", Lambda);
  printEntry(fp,"PRECONDITIONER", precon);
  printEntry(fp,"MULTIGRID SMOOTHER", multigrid_smoother);
  printEntry(fp,"MULTIGRID COARSENING", multigrid_coarsening);
  printEntry(fp,"MULTIGRID CHEBYSHEV DEGREE", multigrid_cheby_degree);
  printEntry(fp,"PARALMOND DEVICE MATRIX TYPE", paralmond_device_matrix_type);
  printEntry(fp,"PARALMOND CYCLE", paralmond_cycle);
  printEntry(fp,"PARALMOND STRENGTH", paralmond_strength);
  printEntry(fp,"PARALMOND RUGESTUBEN STRENGTH THRESHOLD", paralmond_strength_cutoff);
  printEntry(fp,"PARALMOND AGGREGATION", paralmond_aggregation);
  printEntry(fp,"PARALMOND SMOOTHER", paralmond_smoother);
  printEntry(fp,"PARALMOND CHEBYSHEV DEGREE", paralmond_cheby_degree);
  printEntry(fp,"OUTPUT TO FILE", "FALSE");
  printEntry(fp,"VERBOSE", "TRUE");

  fclose(fp);
}

int main(int argc, char **argv){

  for(int degree=1;degree<9;++degree){
    for(double epsy=0.3;epsy<1.1;epsy+=0.7){
      for(int model=1;model<2;++model){
	int maxNx = floor(pow( 6.e6/( pow(degree+1,3)),0.333));
	maxNx = (maxNx>40) ? 40:maxNx;
	maxNx = (maxNx<25) ? 25:maxNx;
	for(int nx=6;nx<maxNx;nx+=6){
	  printSetup("CUDA", degree, nx, nx, nx, epsy);
	  system("../solvers/elliptic/ellipticMain setupHex3D.rc");
	}
      }
    }
  }

}
