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

#include "settings.hpp"

namespace libp {

using std::vector;
using std::string;

setting_t::setting_t(string name_, string val_,
                     string description_, vector<string> options_):
  name{name_},
  val{val_},
  description{description_},
  options{options_} {}

const string& setting_t::getName() const {
  return name;
}

const string& setting_t::getDescription() const {
  return description;
}

const vector<string>& setting_t::getOptions() const {
  return options;
}

void setting_t::updateVal(const string newVal){
  if (!options.size()) {
    val = newVal;
  } else {
    for (size_t i=0;i<options.size();i++) {
      if (newVal==options[i]) {//valid
        val = newVal;
        return;
      }
    }
    stringstream ss;
    ss << "Value: \"" << newVal << "\" "
       << "not valid for setting " << name <<std::endl
       << "Possible values are: { ";
    for (size_t i=0;i<options.size()-1;i++) ss << options[i] << ", ";
    ss << options[options.size()-1] << " }" << std::endl;
    LIBP_FORCE_ABORT(ss.str());
  }
}

bool setting_t::compareVal(const string token) const {
  return !(val.find(token) == std::string::npos);
}

string setting_t::toString() const {
  stringstream ss;

  ss << "Name:  [" << name << "]" << std::endl;
  ss << "Value: " << val << std::endl;

  if (!description.empty())
    ss << "Description: " << description << std::endl;

  if (options.size()) {
    ss << "Possible values: { ";
    for (size_t i=0;i<options.size()-1;i++) ss << options[i] << ", ";
    ss << options[options.size()-1] << " }" << std::endl;
  }

  return ss.str();
}

std::ostream& operator<<(std::ostream& os, const setting_t& setting) {
  os << setting.toString();
  return os;
}

settings_t::settings_t(comm_t _comm):
  comm(_comm) {}

void settings_t::newSetting(const string name, const string val,
                            const string description,
                            const vector<string> options) {
  auto search = settings.find(name);
  if (search == settings.end()) {
    settings[name] = setting_t(name, val, description, options);
    insertOrder.push_back(name);
  } else {
    LIBP_FORCE_ABORT("Setting with name: [" << name << "] already exists.");
  }
}

bool settings_t::hasSetting(const string name) {
  auto search = settings.find(name);
  if (search != settings.end())
    return true;
  else
    return false;
}

void settings_t::changeSetting(const string name, const string newVal) {
  auto search = settings.find(name);
  if (search != settings.end()) {
    setting_t& val = search->second;
    val.updateVal(newVal);
  } else {
    LIBP_FORCE_ABORT("Setting with name: [" << name << "] does not exist.");
  }
}

//input settings from .rc file
void settings_t::readSettingsFromFile(string filename) {
  string line;
  std::ifstream file;

  int rank = comm.rank();

  //only the root rank performs the read
  if (!rank) {
    file.open(filename);
    LIBP_ABORT("Failed to open: " << filename.c_str(),
               !file.is_open());
  }

  string name = "";
  string val  = "";

  //read settings
  int flag;

  if (!rank)
    flag = (getline(file,line)) ? 1 : 0;

  comm.Bcast(flag, 0);

  int MaxLineSize=512;
  memory<char> cline;
  cline.calloc(MaxLineSize+1);

  while (flag) {
    int size;

    if (!rank) {
      size = line.length();
      LIBP_ABORT("Line in settings file is too long: " << line,
                 size>MaxLineSize);
    }
    comm.Bcast(size, 0);

    if (!rank) strcpy(cline.ptr(), line.c_str());
    comm.Bcast(cline, 0, size);

    for(int i=0; i<size; i++){
      char c = cline[i];

      // ignore comments
      if(c == '#') break;

      if(c == '['){ // new setting
        //add current pair if populated
        if (name.length() && val.length()) {
          newSetting(name, val);
          name.clear(); val.clear();
        }

        name=""; val=""; i++;
        while(i < size && cline[i] != ']')
          name += cline[i++];

      // Else add the character
      } else {
        // remove whitespace
        if(isspace(c)) continue;

        val += c;
      }
    }

    if (name.length() && val.length()) {
      newSetting(name, val);
      name.clear(); val.clear();
    }

    if (!rank)
      flag = (getline(file,line)) ? 1 : 0;

    comm.Bcast(flag, 0);
  }

  if (!rank)
    file.close();
}

string settings_t::getSetting(const string name) const {
  auto search = settings.find(name);
  if (search != settings.end()) {
    const setting_t& val = search->second;
    return val.getVal<string>();
  } else {
    LIBP_FORCE_ABORT("Unable to find setting: [" << name << "]");
    return string();
  }
}

bool settings_t::compareSetting(const string name, const string token) const {
  auto search = settings.find(name);
  if (search != settings.end()) {
    const setting_t& val = search->second;
    return val.compareVal(token);
  } else {
    LIBP_FORCE_ABORT("Unable to find setting: [" << name.c_str() << "]");
    return false;
  }
}

void settings_t::report() {
  std::cout << "Settings:\n\n";
  for (size_t i = 0; i < insertOrder.size(); ++i) {
    const string &s = insertOrder[i];
    const setting_t& val = settings[s];
    std::cout << val << std::endl;
  }
}

void settings_t::reportSetting(const string name) const {
  auto search = settings.find(name);
  if (search != settings.end()) {
    const setting_t& val = search->second;
    std::cout << val << std::endl;
  } else {
    LIBP_FORCE_ABORT("Unable to find setting: [" << name.c_str() << "]");
  }
}

} //namespace libp
