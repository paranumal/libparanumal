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

  void print() {
    std::cout << "cheese " << data[0] << " " << data[1] << " " << data[2] << endl;
  }
  
  float squared() const
  {
    float result=0.f;
    for (int i=0; i<3; ++i)
      result+=data[i]*data[i];
    return result;
  }

  //now exclusively for equatorial projection.  Values should already be on the preimage surface.
  v3& normalize()
  {
    float scale = std::sqrt((1 - data[2]*data[2])/(data[0]*data[0]+data[1]*data[1]));
    data[0] = data[0]*scale;
    data[1] = data[1]*scale;
    return *this;
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

using Index=std::uint32_t;

struct Quad
{
  Index vertex[4];
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
  const float Z=0.95f;
  const float X=std::sqrt(1-Z*Z); //autogenerate X from Z
  const float N=0.f;

  static const VertexList vertices=
    {
      {-X,N,Z}, {X,N,Z}, {-X,N,-Z}, {X,N,-Z},
      {N,X,Z}, {N,-X,Z}, {N,X,-Z}, {N,-X,-Z}
    };

  static const QuadList quads=
    {
      {0,2,7,5},{5,7,3,1},{1,3,6,4},{4,6,2,0}
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
  return ((lhs+rhs)/2).normalize();
}

inline ColorPosition split(ColorPosition lhs, ColorPosition rhs)
{
  return{(lhs.color+rhs.color)*0.5f, ((lhs.position+rhs.position)/2).normalize()};
}

template <typename VertexList>
Index vertex_for_edge(Lookup& lookup,
		      VertexList& vertices, Index first, Index second)
{
  Lookup::key_type key(first, second);
  if (key.first>key.second)
    std::swap(key.first, key.second);

  auto inserted=lookup.insert({key, vertices.size()});
  if (inserted.second)
    {
      vertices.push_back(split(vertices[first], vertices[second]));
    }

  return inserted.first->second;
}

template <typename VertexList>
QuadList subdivide_4(VertexList& vertices,
                         QuadList quads)
{
  Lookup lookup;
  QuadList result;

  for (auto&& each:quads)
    {
      std::array<Index, 5> mid; //four on the edge, plus one center
      for (int edge=0; edge<4; ++edge)
	{
	  mid[edge]=vertex_for_edge(lookup, vertices,
				    each.vertex[edge], each.vertex[(edge+1)%4]);
	}

      //midpoint of an arbitrarily chosen edge
      mid[4] = vertex_for_edge(lookup, vertices, mid[0],mid[2]);

      result.push_back({each.vertex[0],mid[0],mid[4],mid[3]});
      result.push_back({each.vertex[1], mid[1],mid[4], mid[0]});
      result.push_back({each.vertex[2], mid[2],mid[4], mid[1]});
      result.push_back({each.vertex[3], mid[3],mid[4], mid[2]});
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

void make_icosphere(int subdivisions, VertexList &vertices, QuadList &quads)
{
  
  vertices=prism::vertices;
  quads=prism::quads;

  for (int i=0; i<subdivisions; ++i)
    {
      quads=subdivide_4(vertices, quads);
    }
}


int main(int argc, char **argv){

  int Nref = (argc==2) ? atoi(argv[1]): 0;
  
  VertexList vertices;
  QuadList quads;
  
  //  IndexedMesh mesh = make_icosphere(2);
  make_icosphere(Nref, vertices, quads);

  // output
  cout << "$MeshFormat" << endl;
  cout << "2.2 0 8" << endl;
  cout << "$EndMeshFormat" << endl;
  cout << "$Nodes" << endl;
  cout <<  vertices.size() << endl;

  cout.precision(15);
  cout << scientific;
  
  // output coordinates of vertices
  
  for(std::vector<int>::size_type i = 0; i != vertices.size(); i++) {
    cout << i+1 << " " 
      <<  vertices[i].data[0] << " "
      <<  vertices[i].data[1] << " "
      <<  vertices[i].data[2] << " "
      <<  endl;
      
  }

  cout << "$EndNodes" << endl;
  cout << "$Elements" << endl;
  cout << quads.size() << endl;

  for(std::vector<int>::size_type i = 0; i != quads.size(); i++) {
    cout
      << i+1 << " "
      << " 3 2 2 6 " 
      <<  quads[i].vertex[0]+1 << " "
      <<  quads[i].vertex[1]+1 << " "
      <<  quads[i].vertex[2]+1 << " "
      <<  quads[i].vertex[3]+1
      <<  endl;
  }
  
  cout << "$EndElements" << endl;

  
  return 0;
}
