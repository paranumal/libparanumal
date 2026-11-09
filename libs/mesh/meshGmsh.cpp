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


#include "mesh/meshGmsh.hpp"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace libp {

namespace gmsh {

namespace {

bool parseLong(char*& p, hlong& value) {
  char* end = nullptr;
  errno = 0;
  const long long v = std::strtoll(p, &end, 10);
  if (end == p || errno == ERANGE) return false;
  value = static_cast<hlong>(v);
  p = end;
  return true;
}

bool parseInt(char*& p, int& value) {
  hlong v = 0;
  if (!parseLong(p, v)) return false;
  value = static_cast<int>(v);
  return true;
}

bool parseDouble(char*& p, double& value) {
  char* end = nullptr;
  const double v = std::strtod(p, &end);
  if (end == p) return false;
  value = v;
  p = end;
  return true;
}

// trimmed first word of a line, e.g. "$Nodes"
std::string firstWord(const char* line) {
  const char* b = line;
  while (*b==' ' || *b=='\t') ++b;
  const char* e = b;
  while (*e && *e!=' ' && *e!='\t' && *e!='\r' && *e!='\n') ++e;
  return std::string(b, e);
}

} //namespace

int reader_t::NumVertices(const int elementType) {
  switch (elementType) {
    case LINE:          return 2;
    case TRIANGLE:      return 3;
    case QUADRILATERAL: return 4;
    case TETRAHEDRON:   return 4;
    case HEXAHEDRON:    return 8;
    default:            return 0;
  }
}

reader_t::reader_t(const std::string _fileName):
  fileName(_fileName) {

  fp = fopen(fileName.c_str(), "r");
  LIBP_ABORT("Cannot open file: " << fileName,
             fp==NULL);

  ReadMeshFormat();

  bool foundNodes = false;
  while (ReadLine()) {
    const std::string section = firstWord(line);
    if (section.size()==0 || section[0]!='$') continue;

    if (section=="$Nodes") {
      if (majorVersion==2) ReadNodesV2();
      else                 ReadNodesV4();
      foundNodes = true;
    } else if (section=="$Elements") {
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " ($Elements found before $Nodes)",
                 !foundNodes);
      fgetpos(fp, &elementsPos);
      foundElements = true;
      break;
    } else if (section=="$Entities" && majorVersion==4) {
      ReadEntities();
    } else if (section=="$PartitionedEntities") {
      LIBP_ABORT("Mesh file " << fileName << " is partitioned by gmsh."
                 " libParanumal partitions meshes itself, please save the mesh"
                 " without partitioning.",
                 true);
    } else {
      SkipSection(section.substr(1));
    }
  }

  LIBP_ABORT("Error reading mesh file: " << fileName
             << " (no $Nodes section)",
             !foundNodes);
  LIBP_ABORT("Error reading mesh file: " << fileName
             << " (no $Elements section)",
             !foundElements);
}

reader_t::~reader_t() {
  if (fp) fclose(fp);
  free(line);
}

bool reader_t::ReadLine() {
  return getline(&line, &lineCapacity, fp) != -1;
}

void reader_t::ReadNextLine() {
  LIBP_ABORT("Error reading mesh file: " << fileName
             << " (unexpected end of file)",
             !ReadLine());
}

void reader_t::SkipSection(const std::string sectionName) {
  const std::string endTag = "$End" + sectionName;
  do {
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (missing " << endTag << ")",
               !ReadLine());
  } while (firstWord(line)!=endTag);
}

void reader_t::ReadMeshFormat() {
  // the first non-empty line must open the $MeshFormat section
  std::string section;
  do {
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (empty file)",
               !ReadLine());
    section = firstWord(line);
  } while (section.size()==0);

  LIBP_ABORT("Error reading mesh file: " << fileName
             << " (expected $MeshFormat, the MSH 1 format is not supported)",
             section!="$MeshFormat");

  ReadNextLine();
  const std::string version = firstWord(line);
  char* p = line + strspn(line, " \t") + version.size();
  int fileType=-1, dataSize=0;
  LIBP_ABORT("Error reading mesh file: " << fileName
             << " (bad $MeshFormat line)",
             !parseInt(p, fileType) || !parseInt(p, dataSize));

  LIBP_ABORT("Mesh file " << fileName << " is a binary .msh file."
             " Only ASCII files are supported, please save it with"
             " Mesh.Binary=0.",
             fileType!=0);

  if (version.compare(0, 2, "2.")==0 || version=="2") {
    majorVersion = 2;
  } else if (version=="4.1") {
    majorVersion = 4;
  } else {
    LIBP_ABORT("Mesh file " << fileName << " uses .msh format version "
               << version << ". Supported versions are 2.2 and 4.1"
               " (e.g. gmsh -format msh41).",
               true);
  }

  SkipSection("MeshFormat");
}

void reader_t::ReadEntities() {
  ReadNextLine();
  char* p = line;
  hlong Nentities[4];
  for (int d=0;d<4;++d) {
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (bad $Entities header)",
               !parseLong(p, Nentities[d]));
  }

  for (int d=0;d<4;++d) {
    // points store x y z, curves/surfaces/volumes store a bounding box
    const int Ncoords = (d==0) ? 3 : 6;
    for (hlong n=0;n<Nentities[d];++n) {
      ReadNextLine();
      p = line;
      int tag=0;
      hlong Nphysical=0, physicalTag=0;
      double x;
      bool ok = parseInt(p, tag);
      for (int c=0;c<Ncoords && ok;++c) ok = parseDouble(p, x);
      ok = ok && parseLong(p, Nphysical);
      if (ok && Nphysical>0) ok = parseLong(p, physicalTag);
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " (bad $Entities entry)",
                 !ok);
      entityPhysicalTag[d][tag] = physicalTag;
    }
  }

  SkipSection("Entities");
}

void reader_t::ReadNodesV2() {
  ReadNextLine();
  char* p = line;
  LIBP_ABORT("Error reading mesh file: " << fileName
             << " (bad $Nodes header)",
             !parseLong(p, Nnodes));

  memory<hlong> tags(Nnodes);
  VX.malloc(Nnodes);
  VY.malloc(Nnodes);
  VZ.malloc(Nnodes);

  for (hlong n=0;n<Nnodes;++n) {
    ReadNextLine();
    p = line;
    double x=0, y=0, z=0;
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (bad node entry)",
               !parseLong(p, tags[n]) || !parseDouble(p, x)
               || !parseDouble(p, y) || !parseDouble(p, z));
    VX[n] = x; VY[n] = y; VZ[n] = z;
  }

  SetNodeOrdering(tags);
  SkipSection("Nodes");
}

void reader_t::ReadNodesV4() {
  ReadNextLine();
  char* p = line;
  hlong Nblocks=0, minTag=0, maxTag=0;
  LIBP_ABORT("Error reading mesh file: " << fileName
             << " (bad $Nodes header)",
             !parseLong(p, Nblocks) || !parseLong(p, Nnodes)
             || !parseLong(p, minTag) || !parseLong(p, maxTag));

  memory<hlong> tags(Nnodes);
  VX.malloc(Nnodes);
  VY.malloc(Nnodes);
  VZ.malloc(Nnodes);

  hlong cnt=0;
  for (hlong b=0;b<Nblocks;++b) {
    ReadNextLine();
    p = line;
    int entityDim=0, entityTag=0, parametric=0;
    hlong NnodesInBlock=0;
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (bad node block header)",
               !parseInt(p, entityDim) || !parseInt(p, entityTag)
               || !parseInt(p, parametric) || !parseLong(p, NnodesInBlock));
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (more nodes than declared in $Nodes header)",
               cnt+NnodesInBlock>Nnodes);

    // node tags of the block, one per line
    for (hlong n=0;n<NnodesInBlock;++n) {
      ReadNextLine();
      p = line;
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " (bad node tag)",
                 !parseLong(p, tags[cnt+n]));
    }
    // then their coordinates (followed by parametric coordinates, if any)
    for (hlong n=0;n<NnodesInBlock;++n) {
      ReadNextLine();
      p = line;
      double x=0, y=0, z=0;
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " (bad node coordinates)",
                 !parseDouble(p, x) || !parseDouble(p, y) || !parseDouble(p, z));
      VX[cnt+n] = x; VY[cnt+n] = y; VZ[cnt+n] = z;
    }
    cnt += NnodesInBlock;
  }

  LIBP_ABORT("Error reading mesh file: " << fileName
             << " (fewer nodes than declared in $Nodes header)",
             cnt!=Nnodes);

  SetNodeOrdering(tags);
  SkipSection("Nodes");
}

/*
  If the node tags are a permutation of 1..Nnodes, store node coordinates
  by tag so a node's index is tag-1 (the v2.2 convention libParanumal has
  always used). Otherwise keep file order and remember tag -> index.
*/
void reader_t::SetNodeOrdering(memory<hlong>& tags) {
  denseNodeTags = true;
  memory<char> seen(Nnodes+1);
  for (hlong n=0;n<=Nnodes;++n) seen[n] = 0;
  for (hlong n=0;n<Nnodes;++n) {
    const hlong tag = tags[n];
    if (tag<1 || tag>Nnodes || seen[tag]) { denseNodeTags = false; break; }
    seen[tag] = 1;
  }

  if (denseNodeTags) {
    bool sorted = true;
    for (hlong n=0;n<Nnodes;++n) {
      if (tags[n]!=n+1) { sorted = false; break; }
    }
    if (!sorted) {
      memory<dfloat> X(Nnodes), Y(Nnodes), Z(Nnodes);
      for (hlong n=0;n<Nnodes;++n) {
        const hlong id = tags[n]-1;
        X[id] = VX[n]; Y[id] = VY[n]; Z[id] = VZ[n];
      }
      VX = X; VY = Y; VZ = Z;
    }
  } else {
    nodeIndex.reserve(Nnodes);
    for (hlong n=0;n<Nnodes;++n) {
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " (duplicate node tag " << tags[n] << ")",
                 !nodeIndex.emplace(tags[n], n).second);
    }
  }
}

hlong reader_t::NodeIndex(const hlong tag) {
  if (denseNodeTags) {
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (element references unknown node " << tag << ")",
               tag<1 || tag>Nnodes);
    return tag-1;
  }
  auto it = nodeIndex.find(tag);
  LIBP_ABORT("Error reading mesh file: " << fileName
             << " (element references unknown node " << tag << ")",
             it==nodeIndex.end());
  return it->second;
}

std::unordered_map<int, hlong> reader_t::CountElements() {
  std::unordered_map<int, hlong> counts;

  fsetpos(fp, &elementsPos);
  ReadNextLine();
  char* p = line;

  if (majorVersion==2) {
    hlong Nelements=0;
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (bad $Elements header)",
               !parseLong(p, Nelements));
    for (hlong n=0;n<Nelements;++n) {
      ReadNextLine();
      p = line;
      hlong tag=0;
      int elementType=0;
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " (bad element entry)",
                 !parseLong(p, tag) || !parseInt(p, elementType));
      ++counts[elementType];
    }
  } else {
    hlong Nblocks=0;
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (bad $Elements header)",
               !parseLong(p, Nblocks));
    for (hlong b=0;b<Nblocks;++b) {
      ReadNextLine();
      p = line;
      int entityDim=0, entityTag=0, elementType=0;
      hlong NelementsInBlock=0;
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " (bad element block header)",
                 !parseInt(p, entityDim) || !parseInt(p, entityTag)
                 || !parseInt(p, elementType) || !parseLong(p, NelementsInBlock));
      counts[elementType] += NelementsInBlock;
      for (hlong n=0;n<NelementsInBlock;++n) ReadNextLine();
    }
  }
  return counts;
}

void reader_t::ForEachElement(const elementCallback_t& f) {
  hlong vertices[8];

  fsetpos(fp, &elementsPos);
  ReadNextLine();
  char* p = line;

  if (majorVersion==2) {
    // elm-number elm-type number-of-tags <tags> node-number-list
    hlong Nelements=0;
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (bad $Elements header)",
               !parseLong(p, Nelements));
    for (hlong n=0;n<Nelements;++n) {
      ReadNextLine();
      p = line;
      hlong tag=0;
      int elementType=0;
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " (bad element entry)",
                 !parseLong(p, tag) || !parseInt(p, elementType));

      const int Nverts = NumVertices(elementType);
      if (Nverts==0) continue;

      // the first tag is the physical group, the number of tags varies
      // (e.g. partitioned meshes add partition tags)
      int Ntags=0;
      hlong physicalTag=0, t=0;
      bool ok = parseInt(p, Ntags) && Ntags>=0;
      for (int i=0;i<Ntags && ok;++i) {
        ok = parseLong(p, t);
        if (i==0) physicalTag = t;
      }
      for (int v=0;v<Nverts && ok;++v) {
        ok = parseLong(p, t);
        if (ok) vertices[v] = NodeIndex(t);
      }
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " (bad element entry)",
                 !ok);
      f(elementType, physicalTag, vertices);
    }
  } else {
    // blocks of elements, each block belongs to one entity
    hlong Nblocks=0;
    LIBP_ABORT("Error reading mesh file: " << fileName
               << " (bad $Elements header)",
               !parseLong(p, Nblocks));
    for (hlong b=0;b<Nblocks;++b) {
      ReadNextLine();
      p = line;
      int entityDim=0, entityTag=0, elementType=0;
      hlong NelementsInBlock=0;
      LIBP_ABORT("Error reading mesh file: " << fileName
                 << " (bad element block header)",
                 !parseInt(p, entityDim) || !parseInt(p, entityTag)
                 || !parseInt(p, elementType) || !parseLong(p, NelementsInBlock)
                 || entityDim<0 || entityDim>3);

      const int Nverts = NumVertices(elementType);
      if (Nverts==0) {
        for (hlong n=0;n<NelementsInBlock;++n) ReadNextLine();
        continue;
      }

      hlong physicalTag=0;
      auto it = entityPhysicalTag[entityDim].find(entityTag);
      if (it!=entityPhysicalTag[entityDim].end()) physicalTag = it->second;

      for (hlong n=0;n<NelementsInBlock;++n) {
        ReadNextLine();
        p = line;
        hlong t=0;
        bool ok = parseLong(p, t); // element tag
        for (int v=0;v<Nverts && ok;++v) {
          ok = parseLong(p, t);
          if (ok) vertices[v] = NodeIndex(t);
        }
        LIBP_ABORT("Error reading mesh file: " << fileName
                   << " (bad element entry)",
                   !ok);
        f(elementType, physicalTag, vertices);
      }
    }
  }
}

} //namespace gmsh

} //namespace libp
