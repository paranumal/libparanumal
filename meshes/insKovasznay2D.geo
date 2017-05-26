r = DefineNumber[ 0.25];
//+
Point(1) = {-0.5, -0.5, 0.0, r};
//+
Point(2) = {-0.5, 1.5,  0.0, r};
//+
Point(3) = {1.0, -0.5,  0.0, r};
//+
Point(4) = {1.0, 1.5, 0.0, r};
//+
Line(1) = {2, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 1};
//+
Line(4) = {1, 2};
//+
Line Loop(5) = {1, 2, 3, 4};
//+
Plane Surface(6) = {5};
//+
Physical Surface("Domain", 9) = {6};
//+
Physical Line("Inflow", 2) = {1, 4, 3};
//+
Physical Line("Outflow", 3) = {2};
//+

