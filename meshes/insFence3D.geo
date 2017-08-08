r1 = DefineNumber[0.050];
r2 = DefineNumber[0.2];
r3 = DefineNumber[0.5];
//+
Point(1) = {0.0, 0.0, 0.0, r1};
Point(2) = {0.0, 0.5, 0.0, r1};
Point(3) = {0.1, 0.5, 0.0, r1};
Point(4) = {0.1, 0.0, 0.0, r1};
//+
Point(5) = {-4.5, 0.0, 0, r2};
Point(6) = {-4.5, 4.5, 0, r3};
Point(7) = {4.5,  0.0, 0, r2};
Point(8) = {4.5,  4.5, 0, r3};
//+
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 5};
Line(5) = {5, 6};
Line(6) = {6, 8};
Line(7) = {8, 7};
Line(8) = {7, 4};
//+
Line Loop(9) = {4, 5, 6, 7, 8, 1, 2, 3};
Plane Surface(10) = {9};
//+
Extrude {0, 0, 2} {Surface{10};}
//+
Physical Surface("Wall",1) = {23, 39, 51, 47, 43};
Physical Surface("Inflow",2) = {27, 52, 31, 10};
Physical Surface("Outflow",3) = {35};
Physical Volume("Domain",9) = {1};
Coherence;
