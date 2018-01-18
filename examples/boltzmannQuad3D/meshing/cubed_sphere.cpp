#include <cmath>
#include <vector>
#include <tuple>
#include <map>
#include <array>

#include <iostream>
using namespace std;

// reference icosasphere meshing from https://github.com/softwareschneiderei/meshing-samples/blob/master/main.cpp

//g++ -g -std=c++11  -o icosasphere icosasphere.cpp  

struct v3
{
  v3(float x, float y, float z)
  {
    data[0]=x;
    data[1]=y;
    data[2]=z;
  }

  v3(float v=0.f) : v3(v,v,v)
  {
  }

  operator const float*() const
  {
    return data;
  }

  v3& operator+=(v3 const& rhs)
  {
    for (int i=0; i<3; ++i)
      data[i]+=rhs[i];
    return *this;
  }
  
  v3& operator-=(v3 const& rhs)
  {
    for (int i=0; i<3; ++i)
      data[i]-=rhs[i];
    return *this;
  }

  v3& operator*=(float rhs)
  {
    for (int i=0; i<3; ++i)
      data[i]*=rhs;
    return *this;
  }

  float squared() const
  {
    float result=0.f;
    for (int i=0; i<3; ++i)
      result+=data[i]*data[i];
    return result;
  }

  float data[3];
};

v3 operator+(v3 lhs, v3 const& rhs)
{
  return lhs+=rhs;
}

v3 operator-(v3 lhs, v3 const& rhs)
{
  return lhs-=rhs;
}

v3 operator*(v3 lhs, float rhs)
{
  return lhs*=rhs;
}

v3 operator/(v3 lhs, float rhs)
{
  return lhs*=(1.f/rhs);
}

v3 normalize(v3 rhs)
{
  return rhs/std::sqrt(rhs.squared());
}

using Index=std::uint32_t;

struct Quad
{
  Index vertex[4];
  bool on_edge[4];
};

struct ColorPosition
{
  v3 color;
  v3 position;
};

using QuadList=std::vector<Quad>;
using VertexList=std::vector<v3>;
using ColorVertexList=std::vector<ColorPosition>;

namespace prism
{
  const float Z=0.5773502691896f;
  const float X=0.5773502691896f;
  const float N=0.5773502691896f;
  
  static const VertexList vertices=
    {
      {-X,N,Z}, {X,N,Z}, {-X,N,-Z}, {X,N,-Z},
      {X,-N,Z}, {-X,-N,Z}, {X,-N,-Z}, {-X,-N,-Z}
    };

  static const QuadList quads=
    {
      {{0,5,7,2},{true,true,true,true}},
      {{5,4,6,7},{true,true,true,true}},
      {{4,1,3,6},{true,true,true,true}},
      {{1,0,2,3},{true,true,true,true}},
      {{0,5,4,1},{true,true,true,true}},
      {{2,7,6,3},{true,true,true,true}}
    };
}

namespace color
{
  v3 red{1.f, 0.f, 0.f};
  v3 green{0.f, 1.f, 0.f};
  v3 blue{0.f, 0.f, 1.f};
  v3 yellow{1.f, 1.f, 0.f};
  v3 cyan{0.f, 1.f, 1.f};
  v3 purple{1.f, 0.f, 1.f};
};


using Lookup=std::map<std::pair<Index, Index>, Index>;

inline v3 split(v3 lhs, v3 rhs)
{
  return normalize(lhs+rhs);
}

inline ColorPosition split(ColorPosition lhs, ColorPosition rhs)
{
  return{(lhs.color+rhs.color)*0.5f, normalize(lhs.position+rhs.position)};
}

template <typename VertexList>
Index vertex_for_edge(Lookup& lookup,
		      VertexList& vertices, std::vector<bool>& mark_verts, Index first, Index second)
{
  Lookup::key_type key(first, second);
  if (key.first>key.second)
    std::swap(key.first, key.second);

  auto inserted=lookup.insert({key, vertices.size()});
  if (inserted.second)
    {
      vertices.push_back(split(vertices[first], vertices[second]));
      mark_verts.push_back(mark_verts[first] && mark_verts[second]);
    }

  return inserted.first->second;
}

template <typename VertexList>
QuadList subdivide_4(VertexList& vertices,
		     std::vector<bool> &mark_verts,
		     QuadList quads)
{
  Lookup lookup;
  QuadList result;

  for (auto&& each:quads)
    {
      std::array<Index, 5> mid; //four on the edge, plus one center
      for (int edge=0; edge<4; ++edge)
	{
	  mid[edge]=vertex_for_edge(lookup, vertices, mark_verts,
				    each.vertex[edge], each.vertex[(edge+1)%4]);
	}

      //midpoint of an arbitrarily chosen edge
      mid[4] = vertex_for_edge(lookup, vertices, mark_verts, mid[0],mid[2]);
      mark_verts.back() = false;

      result.push_back({{each.vertex[0],mid[0],mid[4],mid[3]},{each.on_edge[0],false,false,each.on_edge[3]}});
      result.push_back({{each.vertex[1], mid[1],mid[4], mid[0]},{each.on_edge[1],false,false,each.on_edge[0]}});
      result.push_back({{each.vertex[2], mid[2],mid[4], mid[1]},{each.on_edge[2],false,false,each.on_edge[1]}});
      result.push_back({{each.vertex[3], mid[3],mid[4], mid[2]},{each.on_edge[3],false,false,each.on_edge[2]}});
    }

  return result;
}

int longest_edge(ColorVertexList& vertices, Quad const& quad)
{
  float best=0.f;
  int result=0;
  for (int i=0; i<4; ++i)
    {
      auto a=vertices[quad.vertex[i]].position;
      auto b=vertices[quad.vertex[(i+1)%4]].position;
      float contest =(a-b).squared();
      if (contest>best)
	{
	  best=contest;
	  result=i;
	}
    }
  return result;
}

using IndexedMesh=std::pair<VertexList, QuadList>;
using ColoredIndexMesh=std::pair<ColorVertexList, QuadList>;

void make_icosphere(int subdivisions, VertexList &vertices, std::vector<bool> &mark_verts, QuadList &quads)
{
  
  vertices=prism::vertices;
  quads=prism::quads;
  for (int i = 0; i < vertices.size(); ++i) mark_verts.push_back(true);

  for (int i=0; i<subdivisions; ++i)
    {
      quads=subdivide_4(vertices, mark_verts, quads);
    }
}


int main(int argc, char **argv){

  bool trim_edge = false;

  int Nref = (argc==2) ? atoi(argv[1]): 0;
  
  VertexList vertices;
  QuadList quads;
  std::vector<bool> mark_verts;
  
  //  IndexedMesh mesh = make_icosphere(2);
  make_icosphere(Nref, vertices, mark_verts, quads);

  // output
  cout << "$MeshFormat" << endl;
  cout << "2.2 0 8" << endl;
  cout << "$EndMeshFormat" << endl;
  cout << "$Nodes" << endl;
  if (trim_edge) cout << vertices.size() - 8 - 12 * (pow(2,Nref)-1) << endl;
  else cout <<  vertices.size() << endl;

  cout.precision(15);
  cout << scientific;

  std::vector<int> vertex_offsets;
  
  // output coordinates of vertices
  int j = 0;
  for(std::vector<int>::size_type i = 0; i != vertices.size(); i++) {
    if (!mark_verts[i] || !trim_edge) {
      ++j;
      cout << j << " " 
	   <<  vertices[i].data[0] << " "
	   <<  vertices[i].data[1] << " "
	   <<  vertices[i].data[2] << " "
	   <<  endl;
    }
    vertex_offsets.push_back(j);
  }

  cout.precision(0);
  cout << fixed;
  
  cout << "$EndNodes" << endl;
  cout << "$Elements" << endl;
  if (trim_edge) cout << quads.size() - 12*2*(pow(2,Nref)-1) << endl; 
  else cout << quads.size() << endl;
  
  j = 1;
  for(std::vector<int>::size_type i = 0; i != quads.size(); i++) {
    if (!(quads[i].on_edge[0] || quads[i].on_edge[1] || quads[i].on_edge[2] || quads[i].on_edge[3]) || !trim_edge)
      {
	cout
	  << j << " "
	  << " 3 2 2 6 " 
	  <<  vertex_offsets[quads[i].vertex[0]] << " "
	  <<  vertex_offsets[quads[i].vertex[1]] << " "
	  <<  vertex_offsets[quads[i].vertex[2]] << " "
	  <<  vertex_offsets[quads[i].vertex[3]]
	  <<  endl;
	++j;
      }
  }
  
  cout << "$EndElements" << endl;

  
  return 0;
}
