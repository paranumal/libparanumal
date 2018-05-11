// Base mesh density for outer surfaces
r   = DefineNumber[0.75];
// Factor for inside cylinder faces
fac = DefineNumber[0.1];
// pml Width
pmlWidth = DefineNumber[1.6];
// Define Square Cylinder and Domain boundaries
xmax   = DefineNumber[9.4];
xmin   = DefineNumber[-5.4];
ymax   = DefineNumber[5.4];
ymin   = DefineNumber[-5.4];
zmax   = DefineNumber[5.4];
zmin   = DefineNumber[-5.4];
//
xcmax   = DefineNumber[ 0.5];
xcmin   = DefineNumber[-0.5];
ycmax   = DefineNumber[ 0.5];
ycmin   = DefineNumber[-0.5];
zcmax   = DefineNumber[ 0.5];
zcmin   = DefineNumber[-0.5];
//
Point(1) = {xmin, ymin, zmin, r};
Point(2) = {xmax, ymin, zmin, r};
Point(3) = {xmax, ymax, zmin, r};
Point(4) = {xmin, ymax, zmin, r};
Point(5) = {xmin, ymax, zmin, r};
Point(6) = {xmin, ymin, zmax, r};
Point(7) = {xmax, ymin, zmax, r};
Point(8) = {xmax, ymax, zmax, r};
Point(9) = {xmin, ymax, zmax, r};
Line(1) = {6, 7};
Line(2) = {8, 7};
Line(3) = {8, 9};
Line(4) = {9, 6};
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line(9) = {1, 6};
Line(10) = {4, 9};
Line(11) = {3, 8};
Line(12) = {2, 7};
//
Point(10) = {xcmin, ycmin, zcmin, fac};
Point(11) = {xcmax, ycmin, zcmin, fac};
Point(12) = {xcmax, ycmax, zcmin, fac};
Point(13) = {xcmin, ycmax, zcmin, fac};
Point(14) = {xcmin, ycmin, zcmax, fac};
Point(15) = {xcmax, ycmin, zcmax, fac};
Point(16) = {xcmax, ycmax, zcmax, fac};
Point(17) = {xcmin, ycmax, zcmax, fac};
Line(13) = {14, 15};
Line(14) = {15, 16};
Line(15) = {16, 17};
Line(16) = {17, 14};
Line(17) = {10, 11};
Line(18) = {11, 12};
Line(19) = {12, 13};
Line(20) = {13, 10};
Line(21) = {13, 17};
Line(22) = {12, 16};
Line(23) = {11, 15};
Line(24) = {10, 14};
//
Line Loop(25) = {4, -9, -8, 10};
Plane Surface(26) = {25};
Line Loop(27) = {4, 1, -2, 3};
Plane Surface(28) = {27};
Line Loop(29) = {1, -12, -5, 9};
Plane Surface(30) = {29};
Line Loop(31) = {12, -2, -11, -6};
Plane Surface(32) = {31};
Line Loop(33) = {3, -10, -7, 11};
Plane Surface(34) = {33};
Line Loop(35) = {6, 7, 8, 5};
Plane Surface(36) = {35};
Line Loop(37) = {24, -16, -21, 20};
Plane Surface(38) = {37};
Line Loop(39) = {13, 14, 15, 16};
Plane Surface(40) = {39};
Line Loop(41) = {17, 23, -13, -24};
Plane Surface(42) = {41};
Line Loop(43) = {23, 14, -22, -18};
Plane Surface(44) = {43};
Line Loop(45) = {18, 19, 20, 17};
Plane Surface(46) = {45};
Line Loop(47) = {21, -15, -22, 19};
Plane Surface(48) = {47};
Surface Loop(49) = {40, 42, 46, 44, 48, 38};

 Surface Loop(51) = {34, 28, 26, 30, 32, 36};
 Volume(52) = {51,-49};

Extrude {0, pmlWidth, 0} {Surface{34};}
Extrude {0, 0, -pmlWidth}{Surface{36};}
Extrude {-pmlWidth, 0, 0} {Surface{26};}
Extrude {0, 0, pmlWidth} {Surface{28};}
Extrude {pmlWidth, 0, 0} {Surface{32};}
Extrude {0, -pmlWidth, 0} {Surface{30};}
//
Extrude {0, pmlWidth, 0} {Surface{157, 87, 139, 117};}
Extrude {0, -pmlWidth, 0} {Surface{109, 95, 149, 131};}
Extrude {0, 0, -pmlWidth} {Surface{161, 113};}
Extrude {0, 0, pmlWidth} {Surface{153, 105};}
Extrude {0, -pmlWidth, 0} {Surface{373, 439, 403, 425};}
Extrude {0, pmlWidth, 0} {Surface{381, 395, 417, 447};}

Physical Surface("Wall", 1) = {40, 44, 48, 46, 42, 38};

Physical Surface("Outflow", 3) = {404,294,487, 448, 619, 245, 597, 140, 426, 355, 531, 118, 267, 571, 399, 513, 289, 443, 483, 623, 536, 338, 184, 360, 514, 316, 470, 492, 557, 377, 461, 333, 535, 421, 593, 201, 162, 228, 558, 206, 602, 250, 624, 272, 74, 580, 553, 223, 96, 382, 311, 575, 509, 465};


Physical Volume("Interior",9) = {52};
Physical Volume("XPML",100)   = {55,57};
Physical Volume("YPML",200)   = {53,58};
Physical Volume("ZPML",400)   = {56,54};
Physical Volume("XYPML",300)  = {59,62,63,65};
Physical Volume("XZPML",500)  = {67,68,69,70};
Physical Volume("YZPML",600)  = {60,61,64,66};
Physical Volume("XYZPML",700)  = {75,76,77,78,71,72,73,74};



Coherence;
