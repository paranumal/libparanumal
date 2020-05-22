res = 0.1;
xmin = -1.0; 
xmax =  1.0; 
ymin = -1.0; 
ymax =  1.0; 

Point(1) = {xmin,ymin,0,res};
Point(2) = {xmax,ymin,0,res};
Point(3) = {xmax,ymax,0,res};
Point(4) = {xmin,ymax,0,res};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Physical Surface("Domain", 9) = {6};
Physical Line("Inflow", 1) = {1, 2, 3, 4};


// Point(1) = {xmin,ymin,0,res};
// Point(2) = {xmax,ymin,0,res};
// Point(3) = {xmin,ymax,0,res};
// Line(1) = {1, 2};
// Line(2) = {2, 3};
// Line(3) = {3, 1};
// Line Loop(6) = {1, 2, 3};
// Plane Surface(6) = {6};
// Physical Surface("Domain", 9) = {6};
// Physical Line("Inflow", 1) = {1, 2, 3};