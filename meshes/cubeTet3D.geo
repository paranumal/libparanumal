// Gmsh project created on Mon May 12 11:26:48 2025
SetFactory("OpenCASCADE");
//+
Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};
//+
Box(2) = {-5, -5, -1, 20, 10, 6};

d = 1;
h = d/3;


Field[1] = Box;
//+
Field[1].Thickness = 0.0;
//+
Field[1].VIn = h/5;
//+
Field[1].VOut = 2*h;
//+
Field[1].XMax = 3;
//+
Field[1].XMin = -1;
//+
Field[1].YMax = 1.;
//+
Field[1].YMin = -1;
//+
Field[1].ZMax =  1;
//+
Field[1].ZMin = -1;

Background Field = 1;//+
Surface Loop(3) = {10, 11, 7, 9, 12, 8};
//+
Surface Loop(4) = {2, 3, 5, 1, 4, 6};
//+
Volume(3) = {3, 4};
//+
Physical Volume("Domain", 25) = {2};
//+
Physical Surface("Inflow", 2) = {7, 10, 9};
//+
Physical Surface("Wall", 1) = {2, 4, 1, 3, 5, 11};
//+
Physical Surface("Outflow", 3) = {8};
//+
Physical Surface("Inflow", 2) += {12};
