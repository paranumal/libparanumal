#ifndef __SETUPAIDE
#define __SETUPAIDE

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <sys/types.h>
#include <sys/stat.h>

//#include "headers2d.hpp"
#include "matrix.hpp"

using std::stringstream;
using std::fstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;


class setupAide {
private:
  matrix<string> data;
  matrix<string> keyword;

public:
  setupAide();
  setupAide(string);

  setupAide(const setupAide&);
  setupAide& operator=(const setupAide&);

  string readFile(string);
  void read(string);

  string getArgs(string);

  template <class T>
  int getArgs(string, T&);

  template <class T>
  int getArgs(string, matrix<T>&);

  int getArgs(string, matrix<string>&, string);
};

#include<setupAide.tpp>


#endif

// Combined header/source needs explicit instantiation
