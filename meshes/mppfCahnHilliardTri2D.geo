r = DefineNumber[ 0.125];
xmin = DefineNumber[ 0.0];
xmax = DefineNumber[ 2.0];
ymin = DefineNumber[-1.0];
ymax = DefineNumber[ 1.0];

Point(1) = {xmin, ymin, 0.0, r};
Point(2) = {xmax, ymin, 0.0, r};
Point(3) = {xmax, ymax, 0.0, r};
Point(4) = {xmin, ymax, 0.0, r};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
Line Loop(5) = {1, 2, 3, 4};
//+
Plane Surface(6) = {5};
Physical Surface("Domain", 9) = {6};
Physical Line("Inflow", 2) = {1,2,3,4};


// +Physical Line("Inflow", 2) = {2,3,4};
// +Physical Line("Outflow", 3) = {1};
