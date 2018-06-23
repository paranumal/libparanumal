r0 = DefineNumber[0.025];
r1 = DefineNumber[0.1];
r2 = DefineNumber[0.15];
r3 = DefineNumber[0.15];

//+ Fence Dimensions
fence_max_x = DefineNumber[ 0.05];
fence_min_x = DefineNumber[-0.05];
fence_max_y = DefineNumber[ 0.5];
fence_min_y = DefineNumber[ 0.0];
fence_max_z = DefineNumber[ 0.5];
fence_min_z = DefineNumber[-0.5];

//+ Size of Fine Bounding Box
box_factor_x  = DefineNumber[6.0];
box_factor_y  = DefineNumber[2.0];
box_factor_z  = DefineNumber[2.0];

box_min_x = DefineNumber[box_factor_x*fence_min_x];
box_max_x = DefineNumber[box_factor_x*fence_max_x];
box_min_y = DefineNumber[box_factor_y*fence_min_y];
box_max_y = DefineNumber[box_factor_y*fence_max_y];
box_min_z = DefineNumber[box_factor_z*fence_min_z];
box_max_z = DefineNumber[box_factor_z*fence_max_z];

//+ Size of Wake Volume
wake_max_x = DefineNumber[4.0];
wake_max_y = DefineNumber[1.5*box_max_y];
wake_max_z = DefineNumber[1.5*box_max_z];
wake_min_y = DefineNumber[1.5*box_min_y];
wake_min_z = DefineNumber[1.5*box_min_z];

//+ Size of Domain Without PML
domain_max_x = DefineNumber[4.5];
domain_max_y = DefineNumber[ 2.0];
domain_max_z = DefineNumber[ 2.0];
domain_min_x = DefineNumber[-2.0];
domain_min_y = DefineNumber[ 0.0];
domain_min_z = DefineNumber[-2.0];

// PML Width

pml_width = DefineNumber[0.5];

//+
Point(1) = {fence_min_x, fence_min_y, fence_min_z, r0};
Point(2) = {fence_min_x, fence_max_y, fence_min_z, r0};
Point(3) = {fence_max_x, fence_max_y, fence_min_z, r0};
Point(4) = {fence_max_x, fence_min_y, fence_min_z, r0};

//+
Point(5) = {fence_min_x, fence_min_y, fence_max_z, r0};
Point(6) = {fence_min_x, fence_max_y, fence_max_z, r0};
Point(7) = {fence_max_x, fence_max_y, fence_max_z, r0};
Point(8) = {fence_max_x, fence_min_y, fence_max_z, r0};

//+ Creating Fence Volume 
Line(1) = {5, 8};
Line(2) = {8, 7};
Line(3) = {7, 6};
Line(4) = {6, 5};
Line(5) = {5, 1};
Line(6) = {1, 4};
Line(7) = {4, 8};
Line(8) = {2, 1};
Line(9) = {2, 3};
Line(10) = {3, 4};
Line(11) = {3, 7};
Line(12) = {2, 6};

//+ Creating Bounding Box
Point(9)  = {box_min_x, box_min_y, box_min_z, r0};
Point(10) = {box_min_x, box_max_y, box_min_z, r0};
Point(11) = {box_max_x, box_max_y, box_min_z, r0};
Point(12) = {box_max_x, box_min_y, box_min_z, r0};

Point(13) = {box_min_x, box_min_y, box_max_z, r0};
Point(14) = {box_min_x, box_max_y, box_max_z, r0};
Point(15) = {box_max_x, box_max_y, box_max_z, r0};
Point(16) = {box_max_x, box_min_y, box_max_z, r0};

Line(13) = {13, 16};
Line(14) = {16, 15};
Line(15) = {15, 11};
Line(16) = {11, 10};
Line(17) = {10, 14};
Line(18) = {14, 15};
Line(19) = {14, 13};
Line(20) = {13, 9};
Line(21) = {9, 12};
Line(22) = {12, 16};
Line(23) = {9, 10};
Line(24) = {11, 12};
Line Loop(25) = {20, 21, 22, -13};
Line Loop(26) = {1, -7, -6, -5};
Plane Surface(27) = {25, 26};
Line Loop(28) = {8, -5, -4, -12};
Plane Surface(29) = {28};
Line Loop(30) = {1, 2, 3, 4};
Plane Surface(31) = {30};
Line Loop(32) = {3, -12, 9, 11};
Plane Surface(33) = {32};
Line Loop(34) = {8, 6, -10, -9};
Plane Surface(35) = {34};
Line Loop(36) = {7, 2, -11, 10};
Plane Surface(37) = {36};
Line Loop(38) = {18, -14, -13, -19};
Plane Surface(39) = {38};
Line Loop(40) = {17, 18, 15, 16};
Plane Surface(41) = {40};
Line Loop(42) = {19, 20, 23, 17};
Plane Surface(43) = {42};
Line Loop(44) = {23, -16, 24, -21};
Plane Surface(45) = {44};
Line Loop(46) = {22, 14, 15, 24};
Plane Surface(47) = {46};
Surface Loop(48) = {41, 43, 39, 47, 27, 45, 29, 35, 37, 31, 33};

Volume(1) = {48};

//+ Creating Wake volume
Point(17) = {wake_max_x, wake_min_y, wake_max_z, r1};
Point(18) = {wake_max_x, wake_max_y, wake_max_z, r1};
Point(19) = {wake_max_x, wake_max_y, wake_min_z, r1};
Point(20) = {wake_max_x, wake_min_y, wake_min_z, r1};
Line(49) = {16, 17};
Line(50) = {17, 18};
Line(51) = {18, 15};
Line(52) = {18, 19};
Line(53) = {19, 20};
Line(54) = {20, 17};
Line(55) = {20, 12};
Line(56) = {11, 19};
Line Loop(57) = {49, -54, 55, 22};
Plane Surface(58) = {57};
Line Loop(59) = {54, 50, 52, 53};
Plane Surface(60) = {59};
Line Loop(61) = {51, -14, 49, 50};
Plane Surface(62) = {61};
Line Loop(63) = {56, 53, 55, -24};
Plane Surface(64) = {63};
Line Loop(65) = {15, 56, -52, 51};
Plane Surface(66) = {65};
Surface Loop(67) = {62, 66, 64, 60, 58, 47};
Volume(2) = {67};

//+ Creating domain 
//+ Back Side
Point(21) = {domain_max_x, domain_min_y, domain_max_z, r2};
Point(22) = {domain_max_x, domain_max_y, domain_max_z, r2};
Point(23) = {domain_max_x, domain_min_y, domain_min_z, r2};
Point(24) = {domain_max_x, domain_max_y, domain_min_z, r2};
//+Front Side
Point(25) = {domain_min_x, domain_max_y, domain_max_z, r3};
Point(26) = {domain_min_x, domain_min_y, domain_max_z, r3};
Point(27) = {domain_min_x, domain_max_y, domain_min_z, r3};
Point(28) = {domain_min_x, domain_min_y, domain_min_z, r3};

Line(68) = {26, 21};
Line(69) = {21, 22};
Line(70) = {22, 24};
Line(71) = {24, 23};
Line(72) = {23, 21};
Line(73) = {23, 28};
Line(74) = {28, 26};
Line(75) = {27, 28};
Line(76) = {27, 24};
Line(77) = {27, 25};
Line(78) = {25, 26};
Line(79) = {25, 22};
Line Loop(80) = {68, -72, 73, 74};
Line Loop(81) = {54, -49, -13, 20, 21, -55};
Plane Surface(82) = {80, 81};
Line Loop(83) = {77, 78, -74, -75};
Plane Surface(84) = {83};
Line Loop(85) = {76, 71, 73, -75};
Plane Surface(86) = {85};
Line Loop(87) = {71, 72, 69, 70};
Plane Surface(88) = {87};
Line Loop(89) = {69, -79, 78, 68};
Plane Surface(90) = {89};
Line Loop(91) = {77, 79, 70, -76};
Plane Surface(92) = {91};
Surface Loop(93) = {84, 92, 90, 88, 86, 82, 43, 39, 41, 45, 66, 64, 60, 62};
Volume(3) = {93};


// Now Creating PML with extrusion
// +x pml

Extrude {pml_width, 0, 0} {Surface{88};}
Extrude {-pml_width, 0, 0} {Surface{84};}
//+ 
Extrude {0, 0, -pml_width} {Surface{136, 86, 102};}
Extrude {0, 0, pml_width} {Surface{128, 90, 110};}
//+ 

Extrude {0, pml_width, 0} {Surface{150, 168, 202, 114, 260, 238, 92, 124, 224};}

//+ Domain Volumes
Physical Volume("Interior", 9) = {1, 3, 2};
Physical Volume("XPML",100) = {5, 4};
Physical Volume("YPML",200) = {18};
Physical Volume("ZPML",400) = {10, 7};
Physical Volume("XYPML",300) = {19, 15};
Physical Volume("XZPML",500) = {6, 8, 11, 9};
Physical Volume("YZPML",600) = {13, 17};
Physical Volume("XYZPML",700) = {12, 14, 16, 20};

Physical Surface("Wall",1) = {158, 132, 216, 246, 82, 27, 58, 268, 106, 194, 176, 31, 33, 37, 29, 35};
Physical Surface("Inflow",2) = {220, 137, 154, 282, 440, 466};
Physical Surface("Outflow",3) = {370, 352, 334, 198, 115, 264};
Physical Surface("SlipY",5) = {335, 313, 291, 445, 467, 401, 423, 357, 379};
Physical Surface("SlipZ",6) = {462, 396, 225, 247, 269, 374, 330, 203, 181, 159, 286, 308};

 Coherence;

