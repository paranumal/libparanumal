Mesh.MshFileVersion = 2.2;

d = 1;
h = d/10;

Point(1) = {-3*d,-d,0,h};
Point(2) = {-d,-d,0,h};
Point(3) = {-d,d,0,h};
Point(4) = {-3*d,d,0,h};

Point(5) = {d,-d,0,h};
Point(6) = {d,d,0,h};

Line(7) = {1,2};
Line(8) = {2,3};
Line(9) = {3,4};
Line(10) = {4,1};

Line(11) = {2,5};
Line(12) = {5,6};
Line(13) = {6,3};

Line Loop(14) = {7,8,9,10};
Line Loop(15) = {11,12,13,-8};

Plane Surface(16) = 14;
Plane Surface(17) = 15;

Transfinite Line{7,8,9,10} = d/h;
Transfinite Surface{16};
Recombine Surface{16};//+
Physical Curve("Dirichlet", 1) = {10, 7, 11, 12, 13, 9};
//+
Physical Surface("Domain", 2) = {16, 17};
