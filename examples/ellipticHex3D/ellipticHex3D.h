#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"

void ellipticRunHex3D(mesh3D *mesh);

void ellipticOccaRunHex3D(mesh3D *mesh);

void ellipticSetupHex3D(mesh3D *mesh);

void ellipticVolumeHex3D(mesh3D *mesh);

void ellipticSurfaceHex3D(mesh3D *mesh, dfloat time);

void ellipticUpdateHex3D(mesh3D *mesh, dfloat rka, dfloat rkb);

void ellipticErrorHex3D(mesh3D *mesh, dfloat time);


