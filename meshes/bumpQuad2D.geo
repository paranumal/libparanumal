Mesh.MshFileVersion = 2.2;

SetFactory("OpenCASCADE");

d = 1;
h = d/3;

Point(1) = {-10*d,0,0,h};
Point(2) = {-d/2.,0,0,h};
Point(3) = {-d/2,d,0,h};
Point(4) = { d/2,d,0,h};
Point(5) = { d/2,0,0,h};
Point(6) = {20*d,0,0,h};
Point(7) = {20*d,5*d,0,h};
Point(8) = {-10*d,5*d,0,h};

Line(7) = {1,2};
Line(8) = {2,3};
Line(9) = {3,4};
Line(10) = {4,5};
Line(11) = {5,6};
Line(12) = {6,7};
Line(13) = {7,8};
Line(14) = {8,1};

Line Loop(14) = {7,8,9,10,11,12,13,14};

Plane Surface(15) = 14;

Physical Curve("Wall", 1) = {7,8,9,10,11};
Physical Curve("Inflow", 2) = {14,13};
Physical Curve("Outflow", 3) = {12};
//Physical Curve("Slip", 4) = {13};
//+
Physical Surface("Domain", 20) = {15};

Field[1] = BoundaryLayer;
Field[1].NodesList = {1,2,3,4,5,6,7};
Field[1].CurvesList = {7,8,9,10,11};
Field[1].SizeFar = 0.5*h;
Field[1].Ratio = 2;
Field[1].NbLayers = 8;
Field[1].Thickness = 4*h/2;
Field[1].Quads = 1;
Field[1].FanPointsList = {3,4};
BoundaryLayer Field = 1;

//Mesh.Algorithm = 6;
Recombine Surface {20};
//+
Field[1].NbLayers = 4;
//+
Field[1].Ratio = 4;
//+
Field[1].SizeFar = 0.1;
//+
Field[2] = Min;
//+
Field[2].FieldsList = {1};
//+
Field[1].AnisoMax = 100;
//+
Field[1].Quads = 0;
//+
Field[1].Ratio = 2;
//+
Field[1].NbLayers = 8;
//+
Field[1].IntersectMetrics = 1;
//+
Field[1].Thickness = 0.2;
//+
Field[1].NbLayers = 1;
//+
Field[1].IntersectMetrics = 0;
//+
Field[1].Beta = 2;
//+
Field[1].Size = 0.01;
