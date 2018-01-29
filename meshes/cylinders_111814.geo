lc  = 0.2;
lc1 = 0.4;
lc2 = 2.0;

Point( 1) = { 0.0, -1.0, 0.0 };
Point( 2) = { 0.5, -1.0, 0.0, lc };
Point( 3) = { 0.0, -0.5, 0.0, lc };
Point( 4) = {-0.5, -1.0, 0.0, lc };
Point( 5) = { 0.0, -1.5, 0.0, lc };

Point(11) = { 0.0,  1.0, 0.0 };
Point(12) = { 0.5,  1.0, 0.0, lc };
Point(13) = { 0.0,  1.5, 0.0, lc };
Point(14) = {-0.5,  1.0, 0.0, lc };
Point(15) = { 0.0,  0.5, 0.0, lc };

Point(22) = { 3.0, -1.0, 0.0, lc1 };
Point(32) = { 3.0,  1.0, 0.0, lc1 };

Point(51) = {-20.0, -20.0, 0.0, lc2 };
Point(52) = { 50.0, -20.0, 0.0, lc2 };
Point(53) = { 50.0,  -3.0, 0.0, lc2/3 };
Point(54) = { 50.0,   3.0, 0.0, lc2/3 };
Point(55) = { 50.0,  20.0, 0.0, lc2 };
Point(56) = {-20.0,  20.0, 0.0, lc2 };

Circle( 1) = { 2, 1, 3 };
Circle( 2) = { 3, 1, 4 };
Circle( 3) = { 4, 1, 5 };
Circle( 4) = { 5, 1, 2 };

Circle(11) = {12,11,13 };
Circle(12) = {13,11,14 };
Circle(13) = {14,11,15 };
Circle(14) = {15,11,12 };

Line(21) = {  2, 22 };
Line(22) = { 22, 53 };
Line(23) = { 53, 54 };
Line(24) = { 54, 32 };
Line(25) = { 32, 12 };
Line(26) = { 54, 55 };
Line(27) = { 55, 56 };
Line(28) = { 56, 51 };
Line(29) = { 51, 52 };
Line(30) = { 52, 53 };

Line Loop(1) = { 1, 2, 3, 4, 21, 22, 23, 24, 25,-14,-13,-12,-11,-25,-24,26, 27, 28, 29, 30,-22,-21 };

Plane Surface(1) = { 1 };
Physical Surface("mesh",9) = { 1 };

Physical Line("cyl1",1) = { 1, 2, 3, 4, 11,12,13,14 };
//Physical Line("cyl2",1) = {11,12,13,14 };
Physical Line("inflow",2) = { 28, 29, 27 };
Physical Line("outflow",3) = { 30, 23, 26 };

Field[1] = Attractor;
Field[1].EdgesList = {21,25};

// We then define a Threshold field, which uses the return value of
// the Attractor Field[1] in order to define a simple change in
// element size around the attractors (i.e., around point 5 and line
// 1)
//
// LcMax -                         /------------------
//                               /
//                             /
//                           /
// LcMin -o----------------/
//        |                |       |
//     Attractor       DistMin   DistMax
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 1 / 10;
Field[2].LcMax = 3;
Field[2].DistMin = 1.5;
Field[2].DistMax = 4;

Field[3] = Attractor;
Field[3].EdgesList = {1, 2, 3, 4};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 0.015;
Field[4].LcMax = 2.5;
Field[4].DistMin = 0.1;
Field[4].DistMax = 1.25;

Field[5] = Attractor;
Field[5].EdgesList = {11, 12, 13, 14};

Field[6] = Threshold;
Field[6].IField = 5;
Field[6].LcMin = 0.015;
Field[6].LcMax = 2.5;
Field[6].DistMin = 0.1;
Field[6].DistMax = 1.25;

Field[7] = Min;
Field[7].FieldsList = {2, 4, 6};

Background Field = 7;

Mesh.SecondOrderLinear = 1;
Mesh.ElementOrder = 1;



