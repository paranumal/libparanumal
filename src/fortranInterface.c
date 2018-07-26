#include <string>
#include <cstring>

#include <mpi.h>

#include "setupAide.hpp"
#include "fortranInterface.h"

#define LOWER_BOUND -1

static MPI_Comm *commHandles = NULL;
static int commCount = 0;
static int commMax = 0;
static int commActive = 0;

static setupAide** setupAideHandles = NULL;
static int setupAideCount = 0;
static int setupAideMax = 0;
static int setupAideActive = 0;

int checkHandle(int handle, int low, int up) {
  return (handle > low) && (handle < up);
}

extern "C" {

#define fHolmesInit FORTRAN_NAME(holmesinit,HOLMESINIT)
  void fHolmesInit(int *handle, MPI_Fint *comm, int *err) {
    *err = 1;
    MPI_Comm ccomm = MPI_Comm_f2c(*comm);

    if(commCount == commMax) {
      commMax = commMax / 2 + 1;
      commHandles = (MPI_Comm *) realloc(commHandles, commMax * sizeof(MPI_Comm));
    }

    commHandles[commCount] = ccomm;
    commActive++;
    *handle = commCount++;

    *err = 0;
  }

#define fHolmesFinalize FORTRAN_NAME(holmesfinalize,HOLMESFINALIZE)
  void fHolmesFinalize(int *handle, int *err) {
    *err = 1;

    if(checkHandle(*handle,LOWER_BOUND,commCount)) {
      commHandles[*handle] = 0;
      *handle = LOWER_BOUND;

      commActive--;
      if(commActive == 0) {
        free(commHandles);
        commHandles = NULL;
      }

      *err = 0;
    }
  }

#define fSetupAideCreate FORTRAN_NAME(holmessetupaidecreate,HOLMESSETUPAIDECREATE)
  void fSetupAideCreate(int *handle, int *err) {
    *err = 1;

    if(setupAideCount == setupAideMax) {
      setupAideMax = setupAideMax / 2 + 1;
      setupAideHandles = (setupAide **) realloc(setupAideHandles,
                                              setupAideMax * sizeof(setupAide));
    }

    setupAide *dummy = (setupAide*) calloc(1, sizeof(setupAide));
    setupAideHandles[setupAideCount] = dummy;
    setupAideActive++;
    *handle = setupAideCount++;

    *err = 0;
  }

#define fSetupAideDestroy FORTRAN_NAME(holmessetupaidedestroy,HOLMESSETUPAIDEDESTROY)
  void fSetupAideDestroy(int *handle, int *err) {
    *err = 1;

    if(checkHandle(*handle,LOWER_BOUND,setupAideCount)) {
      free(setupAideHandles[*handle]);
      *handle = LOWER_BOUND;

      setupAideActive--;
      if(setupAideActive == 0) {
        free(setupAideHandles);
      }

      *err = 0;
    }
  }

#define fSetupAideSetArg FORTRAN_NAME(holmessetupaidesetarg,HOLMESSETUPAIDESETARG)
  void fSetupAideSetArg(char *argname, char *argvalue, int *handle, int *err) {
    *err = 1;

    std::string key(argname);
    std::string value(argvalue);

    if(checkHandle(*handle,LOWER_BOUND,setupAideCount)) {
      setupAideHandles[*handle]->setArgs(key, value);
      *err = 0;
    }
  }

#define fSetupAideGetArg FORTRAN_NAME(holmessetupaidegetarg,HOLMESSETUPAIDEGETARG)
  void fSetupAideGetArg(char *argvalue, int *len, char *argname, int *handle, int *err) {
    *err = 1;

    std::string key(argname);
    std::string value;

    if(checkHandle(*handle,LOWER_BOUND,setupAideCount)) {
      setupAideHandles[*handle]->getArgs(key,value);
      std::strcpy(argvalue, value.c_str());
      *len = value.length();
      *err = 0;
    }
  }
}
