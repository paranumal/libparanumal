#include "agmg.h"

void occaTimerTic(occa::device device,std::string name) {
  std::string profilerOn = occa::env::var("OCCA_PROFILE");
  if (profilerOn == "1") {
    device.finish();
    occa::tic(name);
  }
}

void occaTimerToc(occa::device device,std::string name) {
  std::string profilerOn = occa::env::var("OCCA_PROFILE");
  if (profilerOn == "1") {
    device.finish();
    occa::toc(name);
  }
}