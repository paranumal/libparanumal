//#include "headers2d.hpp"

#include "setupAide.hpp"

using std::stringstream;
using std::fstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;

setupAide::setupAide(){}

setupAide::setupAide(string setupFile){
  read(setupFile);
}

setupAide::setupAide(const setupAide& sa){
  *this = sa;
}

setupAide& setupAide::operator = (const setupAide& sa){
  int size = sa.data.size();

  data.resize(size,1);
  keyword.resize(size,1);

  for(int i=1; i<=size; i++){
    data[i]    = sa.data[i];
    keyword[i] = sa.keyword[i];
  }

  return *this;

}

string setupAide::readFile(string filename){
  struct stat statbuf;

  FILE *fh = fopen(filename.c_str(), "r");
  if (fh == 0){
    printf("Failed to open: %s\n", filename.c_str());
    exit(1);
  }

  stat(filename.c_str(), &statbuf);
  char *source = (char *) malloc(statbuf.st_size + 1);
  size_t status = fread(source, statbuf.st_size, 1, fh);
  source[statbuf.st_size] = '\0';

  string ret = source;

  return ret;
}

void setupAide::read(string setupFile){
  vector<string> data2;
  vector<string> keyword2;

  string args = readFile(setupFile);

  int size = args.length();
  string current = "";
  stringstream ss;
  char c;

  for(int i=0; i<size; i++){
    c = args[i];

    // Batch strings together
    if(c == '\'' || c == '"'){
      current += c;
      i++;

      while(i < size && args[i] != c)
	current += args[i++];

      if(i >= size)
	break;

      if( i < (size-1) )
	current += args[i];
    }

    // Batch comments
    else if(c == '/' && i < size && args[i+1] == '*'){
      i += 2;

      while( args[i] != '*' || (i < size && args[i+1] != '/') )
	i++;

      if(i >= size)
	break;

      i++;
    }

    // Removing # comments
    else if(c == '#'){
      i++;

      while(i < size && args[i] != '\n')
	i++;
    }

    // Change \[\] to []
    else if(c == '\\' && i < size && (args[i+1] == '[' || args[i+1] == ']')){
      current += args[i+1];
      i += 2;
    }

    // Split keywords []
    else if(c == '['){
      data2.push_back(current);
      current = "";
      i++;

      while(i < size && args[i] != ']')
	current += args[i++];

      keyword2.push_back(current);
      current = "";
    }

    // Else add the character
    else
      current += c;

    if(i >= (size-1) && current.length())
      data2.push_back(current);

  }

  int argc = (data2.size() - 1);

  data.resize(argc,1);
  keyword.resize(argc,1);

  for(int i=0; i<argc; i++){
    data[i+1]    = data2[i+1];
    keyword[i+1] = keyword2[i];
  }
}

string setupAide::getArgs(string key){
  for(int i=1; i<=keyword.size(); i++)
    if(!( keyword[i].compare(key) ))
      return data[i];

  return "";
}


int setupAide::getArgs(string key, matrix<string>& m, string delimeter){
  string args, current;
  vector<string> argv;
  int argc, size;

  args = getArgs(key);

  size = args.length();

  current = "";

  for(int i=0; i<size; i++){
    while( i < size && delimeter.find(args[i]) == string::npos )
      current += args[i++];

    if(current.length())
      argv.push_back(current);

    current = "";
  }

  argc = argv.size();

  if(!argc){
    printf("Failed to find [%s].\n", key.c_str());
    return 0;
  }

  m.resize(argc,1);

  for(int i=1; i<=argc; i++)
    m[i] = argv[i-1];

  return 1;
}


// Combined header/source needs explicit instantiation
