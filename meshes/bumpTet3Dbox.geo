Mesh.MshFileVersion = 2.2;

SetFactory("OpenCASCADE");

d = 1;
h = d/3;
hz = 5;

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

Field[1] = Box;
//+
Field[1].Thickness = 0.0;
//+
Field[1].VIn = h/5;
//+
Field[1].VOut = h;
//+
Field[1].XMax = 2;
//+
Field[1].XMin = -1;
//+
Field[1].YMax = 1.5;
//+
Field[1].YMin = -2.5;
//+
Field[1].ZMax = hz+1;
//+
Field[1].ZMin = -1;

Field[2] = Box;
//+
Field[2].Thickness = 0;
//+
//Field[2].VIn = h/30;
Field[2].VIn = h/5;
//+
Field[2].VOut = h;
//+
Field[2].XMax = 0.6;
//+
Field[2].XMin = -.6;
//+
Field[2].YMax = 1.1;
//+
Field[2].YMin = -1.1;
//+
Field[2].ZMax = hz+1;
//+
Field[2].ZMin = -1;


Field[3] = Box;
//+
Field[3].Thickness = 0;
//+
//Field[3].VIn = h/30;
Field[3].VIn = h/5;
//+
Field[3].VOut = h;
//+
Field[3].XMax = 20*d;
//+
Field[3].XMin = -10*d;
//+
Field[3].YMax = 0.1;
//+
Field[3].YMin = -0.1;
//+
Field[3].ZMax = hz+1;
//+
Field[3].ZMin = -1;

//+
Field[4] = Min;
//+
Field[4].FieldsList = {1, 2, 3};
//+
Background Field = 4;
//+
Extrude {0, 0, hz} {
  Surface{15}; 
}

//+
Physical Surface("Inflow", 2) = {22, 23, 24, 15};
//+
Physical Surface("Wall", 1) = {17, 18, 20, 16, 19};
//+
Physical Surface("Outflow", 3) = {21};
//+
Physical Volume("Domain", 31) = {1};
