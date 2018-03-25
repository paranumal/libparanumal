//+r1 = DefineNumber[0.0125];
//+r2 = DefineNumber[0.1];
//+r3 = DefineNumber[0.2];
//+r4 = DefineNumber[0.3];

r1 = DefineNumber[0.05];
r2 = DefineNumber[0.1];
r3 = DefineNumber[0.2];
r4 = DefineNumber[0.2];


//+
ffxmin = DefineNumber[-0.05];
ffxmax = DefineNumber[0.05];
ffymin = DefineNumber[0.0];
ffymax = DefineNumber[0.5];

//+
fbfacx = DefineNumber[8.0];
fbfacy = DefineNumber[2.5];

fbxmin = DefineNumber[fbfacx*ffxmin];
fbxmax = DefineNumber[fbfacx*ffxmax];
fbymin = DefineNumber[fbfacy*ffymin];
fbymax = DefineNumber[fbfacy*ffymax];

//+
fdxmin = DefineNumber[-2.5];
fdxmax = DefineNumber[5.0];
fdymin = DefineNumber[0.0];
fdymax = DefineNumber[2.5];

//+
Point(1) = {ffxmin, ffymin, -1.0, r1};
Point(2) = {ffxmin, ffymax, -1.0, r1};
Point(3) = {ffxmax, ffymax, -1.0, r1};
Point(4) = {ffxmax, ffymin, -1.0, r1};

//+
Point(5) = {fbxmin, fbymin, -1.0, r1};
Point(6) = {fbxmin, fbymax, -1.0, r1};
Point(7) = {fbxmax, fbymax, -1.0, r1};
Point(8) = {fbxmax, fbymin, -1.0, r1};

//+
Point(9)  = {fdxmin, fdymin, -1.0, r3};
Point(10) = {fdxmin, fdymax, -1.0, r4};
Point(11) = {fdxmax, fdymax, -1.0, r4};
Point(12) = {fdxmax, fdymin, -1.0, r2};

//+
Point(13) = {fdxmin, fbymax, -1.0, r3};
Point(14) = {fdxmax, fbymax, -1.0, r2};
Point(15) = {fbxmin, fdymax, -1.0, r4};
Point(16) = {fbxmax, fdymax, -1.0, r4};

//+
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 4};
Line(9) = {8, 12};
Line(10) = {12, 14};
Line(11) = {14, 7};
Line(12) = {11, 14};
Line(13) = {16, 7};
Line(14) = {16, 11};
Line(15) = {16, 15};
Line(16) = {15, 6};
Line(17) = {6, 13};
Line(18) = {13, 10};
Line(19) = {10, 15};
Line(20) = {13, 9};
Line(21) = {9, 5};
Line Loop(22) = {5, 6, 7, 8, 1, 2, 3, 4};
Plane Surface(23) = {22};
Line Loop(24) = {7, 9, 10, 11};
Plane Surface(25) = {24};
Line Loop(26) = {11, -13, 14, 12};
Plane Surface(27) = {26};
Line Loop(28) = {16, 6, -13, 15};
Plane Surface(29) = {28};
Line Loop(30) = {16, 17, 18, 19};
Plane Surface(31) = {30};
Line Loop(32) = {17, 20, 21, 5};
Plane Surface(33) = {32};
Recombine Surface {31, 33, 23, 29, 27, 25};
Extrude {0, 0, 2} {
  Surface{31, 29, 27, 25, 23, 33}; Layers{2/r3}; Recombine;
}
Physical Surface("Wall",1) = {154, 158, 150};
Physical Surface("Inflow",2) = {176, 50};
Physical Surface("Outflow",3) = {98, 116};
Physical Surface("YSlip",5) = {54, 76, 94, 112, 146, 162, 180};
Physical Surface("ZSlip",6) = {55, 31, 29, 77, 27, 99, 25, 121, 23, 163, 185, 33};
Physical Volume("Domain",9) = {1, 2, 3, 4, 6, 5};
