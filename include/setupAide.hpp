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

using std::stringstream;
using std::fstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;


class setupAide {
private:
  vector<string> data;
  vector<string> keyword;

public:
  setupAide();
  setupAide(string);

  setupAide(const setupAide&);
  setupAide& operator=(const setupAide&);

  string readFile(string);
  void read(string);

  string getArgs(string);

  void setArgs(string key, string value);

  template <class T>
  int getArgs(string, T&);

  template <class T>
  int getArgs(string, vector<T>&);

  int getArgs(string, vector<string>&, string);


  int compareArgs(string key, string token);

  vector<string> &getData(){ return data; }
  vector<string> &getKeyword() { return keyword; }
};

#include<setupAide.tpp>


#endif

// Combined header/source needs explicit instantiation
