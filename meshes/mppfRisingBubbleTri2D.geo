L = DefineNumber[1.0];
r = DefineNumber[L*0.015];
xmin = DefineNumber[-L/2];
xmax = DefineNumber[ L/2];
ymin = DefineNumber[ 0.0];
ymax = DefineNumber[ 1.5*L];

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
//+Physical Line("Inflow", 2) = {1,2,3,4};
//+Physical Line("Inflow", 2) = {1,3};
Physical Line("YSlip", 5) = {1,3};
Physical Line("XSlip", 4) = {2,4};
