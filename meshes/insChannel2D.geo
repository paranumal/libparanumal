r = DefineNumber[ 0.25];
Point(1) = {0.0, -0.5, 0.0, r};
Point(2) = {5.0, -0.5, 0.0, r};
Point(3) = {5.0,  0.5, 0.0, r};
Point(4) = {0.0,  0.5, 0.0, r};
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
//+
Physical Line("Wall", 1) = {1,3};
Physical Line("Inflow", 2) = {4};
Physical Line("Outflow", 3) = {2};
//+