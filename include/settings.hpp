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

#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <string>
#include <cstring>
#include <vector>
#include <ostream>
#include <sstream>
#include <fstream>
#include "core.hpp"

namespace libp {

class setting_t {
  using string = std::string;
  using stringstream = std::stringstream;

private:
  string name;
  string val;

  string description;
  std::vector<string> options;

public:
  setting_t() = default;
  setting_t(string name_, string val_,
            string description_="", std::vector<string> options_={});

  ~setting_t() = default;

  setting_t(const setting_t& other) = default;
  setting_t& operator=(const setting_t& other) = default;

  const string& getName() const;
  const string& getDescription() const;
  const std::vector<string>& getOptions() const;

  template<typename T>
  T getVal() const {
    stringstream ss(val);
    T arg;
    ss >> arg;
    return arg;
  }

  void updateVal(const string newVal);

  bool compareVal(const string token) const;

  string toString() const;
};

std::ostream& operator<<(std::ostream& os, const setting_t& setting);

class settings_t {
  using string = std::string;
  using stringstream = std::stringstream;

private:
  std::vector<string> insertOrder;

public:
  comm_t comm;
  std::map<string, setting_t> settings;

  settings_t() = default;
  settings_t(comm_t _comm);

  //copy
  settings_t(const settings_t& other)=default;
  settings_t& operator=(const settings_t& other)=default;

  void newSetting(const string name, const string val,
                  const string description="",
                  const std::vector<string> options={});

  bool hasSetting(const string name);

  void changeSetting(const string name, const string newVal);

  //read settings file and update settings
  void readSettingsFromFile(const string filename);

  template <typename T>
  void getSetting(const string name, T& value) const {
    auto search = settings.find(name);
    if (search != settings.end()) {
      const setting_t& val = search->second;
      value = val.getVal<T>();
    } else {
      LIBP_FORCE_ABORT("Unable to find setting: [" << name << "]");
    }
  }

  string getSetting(const string name) const;

  bool compareSetting(const string name, const string token) const;

  void report();

  void reportSetting(const string name) const;
};

} //namespace libp

#endif
