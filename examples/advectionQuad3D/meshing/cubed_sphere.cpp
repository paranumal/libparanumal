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
  v3(double x, double y, double z)
  {
    data[0]=x;
    data[1]=y;
    data[2]=z;
  }

  v3(double v=0.f) : v3(v,v,v)
  {
  }

  operator const double*() const
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

  v3& operator*=(double rhs)
  {
    for (int i=0; i<3; ++i)
      data[i]*=rhs;
    return *this;
  }

  double squared() const
  {
    double result=0.f;
    for (int i=0; i<3; ++i)
      result+=data[i]*data[i];
    return result;
  }

  double data[3];
};

v3 operator+(v3 lhs, v3 const& rhs)
{
  return lhs+=rhs;
}

v3 operator-(v3 lhs, v3 const& rhs)
{
  return lhs-=rhs;
}

v3 operator*(v3 lhs, double rhs)
{
  return lhs*=rhs;
}

v3 operator/(v3 lhs, double rhs)
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
  Index vertex[4]; //the physical location of the quad
  bool on_edge[4]; //marks the edge
  int face; //which face of the cube the quad is on
  int dist[4]; //distance to edge of cube
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
  const double Z=0.5773502691896f;
  const double X=0.5773502691896f;
  const double N=0.5773502691896f;
  
  static const VertexList vertices=
    {
      {-X,N,Z}, {X,N,Z}, {-X,N,-Z}, {X,N,-Z},
      {X,-N,Z}, {-X,-N,Z}, {X,-N,-Z}, {-X,-N,-Z}
    };

  static const QuadList quads=
    {
      {{2,7,5,0},{true,true,true,true},3,{0,0,0,0}},
      {{7,6,4,5},{true,true,true,true},4,{0,0,0,0}},
      {{6,3,1,4},{true,true,true,true},1,{0,0,0,0}},
      {{3,2,0,1},{true,true,true,true},2,{0,0,0,0}},
      {{4,1,0,5},{true,true,true,true},5,{0,0,0,0}},
      {{7,2,3,6},{true,true,true,true},6,{0,0,0,0}}
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
  //return normalize(lhs+rhs);
  return (lhs+rhs)/2;
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

      if (each.dist[0] == 0 && each.dist[1] == 0 && each.dist[2] == 0 && each.dist[3] == 0) {
	result.push_back({{each.vertex[0],mid[0],mid[4],mid[3]},{each.on_edge[0],false,false,each.on_edge[3]},each.face,{0,0,1,0}});
	result.push_back({{mid[0],each.vertex[1], mid[1],mid[4]},{each.on_edge[1],false,false,each.on_edge[0]},each.face,{0,0,0,1}});
	result.push_back({{mid[4], mid[1],each.vertex[2], mid[2]},{each.on_edge[2],false,false,each.on_edge[1]},each.face,{1,0,0,0}});
	result.push_back({{ mid[3],mid[4], mid[2],each.vertex[3]},{each.on_edge[3],false,false,each.on_edge[2]},each.face,{0,1,0,0}});
      }
      else {
	int max4 = max(each.dist[0]+each.dist[2],each.dist[1]+each.dist[3]);
	
	result.push_back({{each.vertex[0],mid[0],mid[4],mid[3]},{each.on_edge[0],false,false,each.on_edge[3]},each.face,{2*each.dist[0],each.dist[0] + each.dist[1],max4,each.dist[0]+each.dist[3]}});
	result.push_back({{mid[0],each.vertex[1], mid[1],mid[4]},{each.on_edge[1],false,false,each.on_edge[0]},each.face,{each.dist[0] + each.dist[1],2*each.dist[1],each.dist[1]+each.dist[2],max4}});
	result.push_back({{mid[4], mid[1],each.vertex[2], mid[2]},{each.on_edge[2],false,false,each.on_edge[1]},each.face,{max4,each.dist[1]+each.dist[2],2*each.dist[2],each.dist[2]+each.dist[3]}});
	result.push_back({{ mid[3],mid[4], mid[2],each.vertex[3]},{each.on_edge[3],false,false,each.on_edge[2]},each.face,{each.dist[3]+each.dist[0],max4,each.dist[2]+each.dist[3],2*each.dist[3]}});
      }
    }

  return result;
}

int longest_edge(ColorVertexList& vertices, Quad const& quad)
{
  double best=0.f;
  int result=0;
  for (int i=0; i<4; ++i)
    {
      auto a=vertices[quad.vertex[i]].position;
      auto b=vertices[quad.vertex[(i+1)%4]].position;
      double contest =(a-b).squared();
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

  bool trim_edge = false; //set to true to remove cube edge entirely
  bool mark_edge = true; //adds a tag to mark cube face matching each element face

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
    //
    if (!trim_edge || !(quads[i].on_edge[0] ||
			quads[i].on_edge[1] ||
			quads[i].on_edge[2] ||
			quads[i].on_edge[3]))
      {

	if (mark_edge) {
	  int min_dist = 1e9;
	  for (int j = 0; j < 4; ++j) {
	    min_dist = min(min_dist,quads[i].dist[j]);
	  }
	  cout
	    << j << " "
	    << " 3 4 2 6 "
	    <<  quads[i].face << " "
	    <<  min_dist << " "
	    <<  vertex_offsets[quads[i].vertex[0]] << " "
	    <<  vertex_offsets[quads[i].vertex[1]] << " "
	    <<  vertex_offsets[quads[i].vertex[2]] << " "
	    <<  vertex_offsets[quads[i].vertex[3]]
	    <<  endl;
	     }
	else
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
