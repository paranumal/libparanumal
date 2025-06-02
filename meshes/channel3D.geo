// Gmsh project created on Mon Jun  2 09:27:48 2025
SetFactory("OpenCASCADE");
//+
Box(1) = {-1, -1, -1, 8, 2, 2};
//+
Surface Loop(2) = {5, 1, 3, 6, 2, 4};
//+
Volume(2) = {2};
//+
Physical Surface("Inflow", 2) = {5, 3, 4, 1, 6};
//+
Physical Surface("Outflow", 3) = {2};
//+
Physical Volume("Domain", 10) = {1};
