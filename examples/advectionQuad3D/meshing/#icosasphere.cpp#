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

struct Triangle
{
  Index vertex[3];
};

struct ColorPosition
{
  v3 color;
  v3 position;
};

using TriangleList=std::vector<Triangle>;
using VertexList=std::vector<v3>;
using ColorVertexList=std::vector<ColorPosition>;

namespace icosahedron
{
  const float X=.525731112119133606f;
  const float Z=.850650808352039932f;
  const float N=0.f;

  static const VertexList vertices=
    {
      {-X,N,Z}, {X,N,Z}, {-X,N,-Z}, {X,N,-Z},
      {N,Z,X}, {N,Z,-X}, {N,-Z,X}, {N,-Z,-X},
      {Z,X,N}, {-Z,X, N}, {Z,-X,N}, {-Z,-X, N}
    };

  static const TriangleList triangles=
    {
      {0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
      {8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
      {7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
      {6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
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
TriangleList subdivide_4(VertexList& vertices,
                         TriangleList triangles)
{
  Lookup lookup;
  TriangleList result;

  for (auto&& each:triangles)
    {
      std::array<Index, 3> mid;
      for (int edge=0; edge<3; ++edge)
	{
	  mid[edge]=vertex_for_edge(lookup, vertices,
				    each.vertex[edge], each.vertex[(edge+1)%3]);
	}

      result.push_back({each.vertex[0], mid[0], mid[2]});
      result.push_back({each.vertex[1], mid[1], mid[0]});
      result.push_back({each.vertex[2], mid[2], mid[1]});
      result.push_back({mid[0], mid[1], mid[2]});
    }

  return result;
}

int longest_edge(ColorVertexList& vertices, Triangle const& triangle)
{
  float best=0.f;
  int result=0;
  for (int i=0; i<3; ++i)
    {
      auto a=vertices[triangle.vertex[i]].position;
      auto b=vertices[triangle.vertex[(i+1)%3]].position;
      float contest =(a-b).squared();
      if (contest>best)
	{
	  best=contest;
	  result=i;
	}
    }
  return result;
}

TriangleList subdivide_2(ColorVertexList& vertices,
                         TriangleList triangles)
{
  Lookup lookup;
  TriangleList result;

  for (auto&& each:triangles)
    {
      auto edge=longest_edge(vertices, each);
      Index mid=vertex_for_edge(lookup, vertices,
				each.vertex[edge], each.vertex[(edge+1)%3]);

      result.push_back({each.vertex[edge],
	    mid, each.vertex[(edge+2)%3]});

      result.push_back({each.vertex[(edge+2)%3],
	    mid, each.vertex[(edge+1)%3]});
    }

  return result;
}

using IndexedMesh=std::pair<VertexList, TriangleList>;
using ColoredIndexMesh=std::pair<ColorVertexList, TriangleList>;

void make_icosphere(int subdivisions, VertexList &vertices, TriangleList &triangles)
{
  //  VertexList vertices=icosahedron::vertices;
  //  TriangleList triangles=icosahedron::triangles;

  vertices=icosahedron::vertices;
  triangles=icosahedron::triangles;

  for (int i=0; i<subdivisions; ++i)
    {
      triangles=subdivide_4(vertices, triangles);
    }

  //  return{vertices, triangles};
}


int main(int argc, char **argv){

  int Nref = (argc==2) ? atoi(argv[1]): 0;
  
  VertexList vertices;
  TriangleList triangles;
  
  //  IndexedMesh mesh = make_icosphere(2);
  make_icosphere(Nref, vertices, triangles);

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
      <<  vertices[i].data[2] << endl;
      
  }

  cout << "$EndNodes" << endl;
  cout << "$Elements" << endl;
  cout << triangles.size() << endl;

  for(std::vector<int>::size_type i = 0; i != triangles.size(); i++) {
    cout
      << i+1 << " "
      << " 2 2 0 1 " 
      <<  triangles[i].vertex[0]+1 << " "
      <<  triangles[i].vertex[2]+1 << " "
      <<  triangles[i].vertex[1]+1 << " "
      <<  endl;
  }
  
  cout << "$EndElements" << endl;

  
  return 0;
}
