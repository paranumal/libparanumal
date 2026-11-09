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


#ifndef MESH_GMSH_HPP
#define MESH_GMSH_HPP 1

#include "core.hpp"
#include <functional>
#include <unordered_map>

namespace libp {

namespace gmsh {

// gmsh element type codes used by libParanumal
constexpr int LINE          = 1;
constexpr int TRIANGLE      = 2;
constexpr int QUADRILATERAL = 3;
constexpr int TETRAHEDRON   = 4;
constexpr int HEXAHEDRON    = 5;

/*
  Reader for ASCII gmsh .msh files in format version 2.2 or 4.1.

  The constructor reads the node coordinates. Elements are then streamed
  from the file with CountElements and ForEachElement, so callers only have
  to store the elements they keep.

  Node indices handed to callers are 0-based. When the node tags in the
  file are exactly 1..Nnodes (what gmsh writes by default) the index of a
  node is its tag minus one, as in the previous v2.2-only reader. Otherwise
  nodes are numbered in the order they appear in the file.

  The physical tag of an element is the first tag in its tag list (v2.2)
  or the first physical tag of the entity its block belongs to (v4.1). It is
  0 when the element has no physical group.
*/
class reader_t {
public:
  using elementCallback_t = std::function<void(const int elementType,
                                               const hlong physicalTag,
                                               const hlong* vertices)>;

  reader_t(const std::string fileName);
  ~reader_t();

  reader_t(const reader_t&) = delete;
  reader_t& operator=(const reader_t&) = delete;

  int MajorVersion() const { return majorVersion; }

  hlong Nnodes=0;
  memory<dfloat> VX, VY, VZ;

  // number of elements of each gmsh element type in the file
  std::unordered_map<int, hlong> CountElements();

  // call f(elementType, physicalTag, vertices) for every line, triangle,
  // quadrilateral, tetrahedron and hexahedron in the file, in file order.
  // Other element types (points, high order elements) are skipped.
  void ForEachElement(const elementCallback_t& f);

  // number of vertices of a supported gmsh element type, 0 otherwise
  static int NumVertices(const int elementType);

private:
  std::string fileName;
  FILE *fp=nullptr;
  char *line=nullptr;
  size_t lineCapacity=0;

  int majorVersion=0;
  fpos_t elementsPos;
  bool foundElements=false;

  // v4.1: first physical tag of each (dim, entity tag)
  std::unordered_map<int, hlong> entityPhysicalTag[4];

  // node tag -> node index, used only when tags are not 1..Nnodes
  std::unordered_map<hlong, hlong> nodeIndex;
  bool denseNodeTags=true;

  bool ReadLine();
  void ReadNextLine();
  void SkipSection(const std::string sectionName);
  void ReadMeshFormat();
  void ReadEntities();
  void ReadNodesV2();
  void ReadNodesV4();
  void SetNodeOrdering(memory<hlong>& tags);
  hlong NodeIndex(const hlong tag);
};

} //namespace gmsh

} //namespace libp

#endif
