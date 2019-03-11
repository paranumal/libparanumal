#include <cmath>
#include <vector>
#include <tuple>
#include <map>
#include <array>
#include <algorithm>

#include <iostream>
using namespace std;

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
    return rhs;///std::sqrt(rhs.squared());
}

struct Quad
{
  int vertex[4]; //the physical location of the quad
  int face; //which face of the cube the quad is on
  int dist; //distance to edge of cube
};

using QuadList=std::vector<Quad>;
using VertexList=std::vector<v3>;
using CoordinateList=std::vector<double>;
using dedup=std::map<std::array<int,3>, int>;

void make_coordinates(int subdivisions,CoordinateList *axis) {
  double start = -1/sqrt(3);
  double end = 1/sqrt(3);
  double step = (end - start)/subdivisions;
  for (int i = 0; i <= subdivisions; ++i)
    axis->push_back(start + i*step);
  return;
}

void  make_grid(int subdivisions, VertexList &vertices, QuadList& quads) {
  CoordinateList axis;
  dedup indices;
  double face_offset = 1/sqrt(3);
  make_coordinates(subdivisions,&axis);
  for (int xi = 0; xi < subdivisions; ++xi) {
    for (int yi = 0; yi < subdivisions; ++yi) {
      //face 1
      v3 e1 = {face_offset,axis[xi],axis[yi]};
      v3 e2 = {face_offset,axis[xi+1],axis[yi]};
      v3 e3 = {face_offset,axis[xi+1],axis[yi+1]};
      v3 e4 = {face_offset,axis[xi],axis[yi+1]};

      std::array<int,3> eLat1 = {subdivisions,xi,yi};
      std::array<int,3> eLat2 = {subdivisions,xi+1,yi};
      std::array<int,3> eLat3 = {subdivisions,xi+1,yi+1};
      std::array<int,3> eLat4 = {subdivisions,xi,yi+1};

      auto map1 = indices.insert({eLat1,vertices.size()});
      if (map1.second) vertices.push_back(e1);
      auto map2 = indices.insert({eLat2,vertices.size()});
      if (map2.second) vertices.push_back(e2);
      auto map3 = indices.insert({eLat3,vertices.size()});
      if (map3.second) vertices.push_back(e3);
      auto map4 = indices.insert({eLat4,vertices.size()});
      if (map4.second) vertices.push_back(e4);
      
      int faceNumber = 1;

      int dists1[] = {xi,yi,subdivisions - 1 - xi,subdivisions - 1 - yi};
      int min_dist = *min_element(dists1,dists1+4);
      
      quads.push_back({{map1.first->second,map2.first->second,map3.first->second,map4.first->second},faceNumber,min_dist});

      //face 2
      /*      e1 = {axis[subdivisions-xi],face_offset,axis[yi]};
      e2 = {axis[subdivisions-xi-1],face_offset,axis[yi]};
      e3 = {axis[subdivisions-xi-1],face_offset,axis[yi+1]};
      e4 = {axis[subdivisions-xi],face_offset,axis[yi+1]};

      eLat1 = {subdivisions-xi,subdivisions,yi};
      eLat2 = {subdivisions-xi-1,subdivisions,yi};
      eLat3 = {subdivisions-xi-1,subdivisions,yi+1};
      eLat4 = {subdivisions-xi,subdivisions,yi+1};
      
      map1 = indices.insert({eLat1,vertices.size()});
      if (map1.second) vertices.push_back(e1);
      map2 = indices.insert({eLat2,vertices.size()});
      if (map2.second) vertices.push_back(e2);
      map3 = indices.insert({eLat3,vertices.size()});
      if (map3.second) vertices.push_back(e3);
      map4 = indices.insert({eLat4,vertices.size()});
      if (map4.second) vertices.push_back(e4);
      
      faceNumber = 2;

      int dists2[] = {xi,yi,subdivisions - 1 - xi,subdivisions - 1 - yi};
      min_dist = *min_element(dists2,dists2+4);

      quads.push_back({{map1.first->second,map2.first->second,map3.first->second,map4.first->second},faceNumber,min_dist});

      //face 3
      e1 = {-1*face_offset,axis[subdivisions-xi],axis[yi]};
      e2 = {-1*face_offset,axis[subdivisions-xi-1],axis[yi]};
      e3 = {-1*face_offset,axis[subdivisions-xi-1],axis[yi+1]};
      e4 = {-1*face_offset,axis[subdivisions-xi],axis[yi+1]};

      eLat1 = {0,subdivisions-xi,yi};
      eLat2 = {0,subdivisions-xi-1,yi};
      eLat3 = {0,subdivisions-xi-1,yi+1};
      eLat4 = {0,subdivisions-xi,yi+1};
      
      map1 = indices.insert({eLat1,vertices.size()});
      if (map1.second) vertices.push_back(e1);
      map2 = indices.insert({eLat2,vertices.size()});
      if (map2.second) vertices.push_back(e2);
      map3 = indices.insert({eLat3,vertices.size()});
      if (map3.second) vertices.push_back(e3);
      map4 = indices.insert({eLat4,vertices.size()});
      if (map4.second) vertices.push_back(e4);
      
      faceNumber = 3;

      int dists3[] = {xi,yi,subdivisions - 1 - xi,subdivisions - 1 - yi};
      min_dist = *min_element(dists3,dists3+4);
      
      quads.push_back({{map1.first->second,map2.first->second,map3.first->second,map4.first->second},faceNumber,min_dist});

      //face 4
      e1 = {axis[xi],-1*face_offset,axis[yi]};
      e2 = {axis[xi+1],-1*face_offset,axis[yi]};
      e3 = {axis[xi+1],-1*face_offset,axis[yi+1]};
      e4 = {axis[xi],-1*face_offset,axis[yi+1]};

      eLat1 = {xi,0,yi};
      eLat2 = {xi+1,0,yi};
      eLat3 = {xi+1,0,yi+1};
      eLat4 = {xi,0,yi+1};
      
      map1 = indices.insert({eLat1,vertices.size()});
      if (map1.second) vertices.push_back(e1);
      map2 = indices.insert({eLat2,vertices.size()});
      if (map2.second) vertices.push_back(e2);
      map3 = indices.insert({eLat3,vertices.size()});
      if (map3.second) vertices.push_back(e3);
      map4 = indices.insert({eLat4,vertices.size()});
      if (map4.second) vertices.push_back(e4);
      
      faceNumber = 4;

      int dists4[] = {xi,yi,subdivisions - 1 - xi,subdivisions - 1 - yi};
      min_dist = *min_element(dists4,dists4+4);
      
      quads.push_back({{map1.first->second,map2.first->second,map3.first->second,map4.first->second},faceNumber,min_dist});

      //face 5
      e1 = {axis[subdivisions-xi],axis[yi],face_offset};
      e2 = {axis[subdivisions-xi],axis[yi+1],face_offset};
      e3 = {axis[subdivisions-xi-1],axis[yi+1],face_offset};
      e4 = {axis[subdivisions-xi-1],axis[yi],face_offset};

      eLat1 = {subdivisions-xi,yi,subdivisions};
      eLat2 = {subdivisions-xi,yi+1,subdivisions};
      eLat3 = {subdivisions-xi-1,yi+1,subdivisions};
      eLat4 = {subdivisions-xi-1,yi,subdivisions};
      
      map1 = indices.insert({eLat1,vertices.size()});
      if (map1.second) vertices.push_back(e1);
      map2 = indices.insert({eLat2,vertices.size()});
      if (map2.second) vertices.push_back(e2);
      map3 = indices.insert({eLat3,vertices.size()});
      if (map3.second) vertices.push_back(e3);
      map4 = indices.insert({eLat4,vertices.size()});
      if (map4.second) vertices.push_back(e4);
      
      faceNumber = 5;

      int dists5[] = {xi,yi,subdivisions - 1 - xi,subdivisions - 1 - yi};
      min_dist = *min_element(dists5,dists5+4);
      
      quads.push_back({{map1.first->second,map2.first->second,map3.first->second,map4.first->second},faceNumber,min_dist});

      //face 6
      e1 = {axis[xi],axis[yi],-1*face_offset};
      e2 = {axis[xi],axis[yi+1],-1*face_offset};
      e3 = {axis[xi+1],axis[yi+1],-1*face_offset};
      e4 = {axis[xi+1],axis[yi],-1*face_offset};

      eLat1 = {xi,yi,0};
      eLat2 = {xi,yi+1,0};
      eLat3 = {xi+1,yi+1,0};
      eLat4 = {xi+1,yi,0};
      
      map1 = indices.insert({eLat1,vertices.size()});
      if (map1.second) vertices.push_back(e1);
      map2 = indices.insert({eLat2,vertices.size()});
      if (map2.second) vertices.push_back(e2);
      map3 = indices.insert({eLat3,vertices.size()});
      if (map3.second) vertices.push_back(e3);
      map4 = indices.insert({eLat4,vertices.size()});
      if (map4.second) vertices.push_back(e4);

      faceNumber = 6;

      int dists6[] = {xi,yi,subdivisions - 1 - xi,subdivisions - 1 - yi};
      min_dist = *min_element(dists6,dists6+4);

      quads.push_back({{map1.first->second,map2.first->second,map3.first->second,map4.first->second},faceNumber,min_dist});*/
    }
  }
}

int main(int argc, char **argv){

  int Nref = (argc==2) ? atoi(argv[1]): 0;
  
  VertexList vertices;
  QuadList quads;
  
  make_grid(Nref, vertices, quads);

  // output
  cout << "$MeshFormat" << endl;
  cout << "2.2 0 8" << endl;
  cout << "$EndMeshFormat" << endl;
  cout << "$Nodes" << endl;
  cout <<  vertices.size() << endl;

  cout.precision(15);
  cout << scientific;

  std::vector<int> vertex_offsets;
  
  // output coordinates of vertices
  for(int i = 0; i < vertices.size(); i++) {
    cout << i+1 << " " 
	 <<  normalize(vertices[i]).data[0] << " "
	 <<  normalize(vertices[i]).data[1] << " "
	 <<  normalize(vertices[i]).data[2] << " "
	 <<  endl;
  }

  cout.precision(0);
  cout << fixed;
  
  cout << "$EndNodes" << endl;
  cout << "$Elements" << endl;
  cout << quads.size() << endl;
  
  for(int i = 0; i < quads.size(); i++) {
    cout
      << i+1 << " "
      << " 3 4 2 6 "
      <<  quads[i].face << " "
      <<  quads[i].dist << " "
      <<  quads[i].vertex[0]+1 << " "
      <<  quads[i].vertex[1]+1 << " "
      <<  quads[i].vertex[2]+1 << " "
      <<  quads[i].vertex[3]+1
      <<  endl;
  }
  
  cout << "$EndElements" << endl;

  
  return 0;
}
