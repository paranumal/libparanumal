cl__1 = 1.0;
xmax = DefineNumber[1.0];
xmin = DefineNumber[0.0];

// xn = DefineNumber[3];
// yn = DefineNumber[3];
// zn = DefineNumber[3];

Point(1) = {xmin, xmin, xmin, cl__1};
Point(2) = {xmax, xmin, xmin, cl__1};
Point(3) = {xmax, xmax, xmin, cl__1};
Point(4) = {xmin, xmax, xmin, cl__1};
Point(5) = {xmin, xmin, xmax, cl__1};
Point(6) = {xmax, xmin, xmax, cl__1};
Point(7) = {xmax, xmax, xmax, cl__1};
Point(8) = {xmin, xmax, xmax, cl__1};
//
Line(9)  = {1, 2};
Line(10) = {2, 3};
Line(11) = {3, 4};
Line(12) = {4, 1};
Line(13) = {5, 6};
Line(14) = {6, 7};
Line(15) = {7, 8};
Line(16) = {8, 5};
Line(17) = {1, 5};
Line(18) = {2, 6};
Line(19) = {3, 7};
Line(20) = {4, 8};
//
Transfinite Line {9, 13, 11, 15}  = (2^xn+1) Using Progression 1;
Transfinite Line {12, 10, 14, 16} = (2^yn+1) Using Progression 1;
Transfinite Line {17, 18, 20, 19} = (2^zn+1) Using Progression 1;
//
Line Loop(22) = {16, 13, 14, 15};
Plane Surface(22) = {22};
Transfinite Surface {22} = {6,7,8,5};
Recombine Surface {22};
Line Loop(24) = {14, -19, -10, 18};
Plane Surface(24) = {24};
Transfinite Surface {24} = {6,2,3,7};
Recombine Surface {24};
Line Loop(26) = {10, 11, 12, 9};
Plane Surface(26) = {26};
Transfinite Surface {26};
Recombine Surface {26};
Line Loop(28) = {12, 17, -16, -20};
Plane Surface(28) = {28};
Transfinite Surface {28} = {5,8,4,1};
Recombine Surface {28};
Line Loop(30) = {20, -15, -19, 11};
Plane Surface(30) = {30};
Transfinite Surface {30};
Recombine Surface {30};
Line Loop(32) = {18, -13, -17, 9};
Plane Surface(32) = {32};
Transfinite Surface {32};
Recombine Surface {32};
Surface Loop(34) = {22, 28, 26, 24, 30, 32};
Volume(34) = {34};
Transfinite Volume {34} = {5,6,2,1,8,7,3,4};
//+Physical Surface("Inflow") = {22, 24, 26, 28, 30, 32};
//+Physical Volume("Domain") = {34};

