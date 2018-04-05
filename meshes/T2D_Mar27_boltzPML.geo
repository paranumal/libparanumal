
lc = 1.00;

// x-bounds
xmin  = -5.0;
xmax  =  13.0;
pmlx1 =  2.0;
pmlx2 =  2.0;

// y-bounds
ymin  =  -6.0;
ymax  =  6.0;
pmly1 =  2.0;
pmly2 =  2.0;

// z-bounds
zall = 0.0;

h = 1.0;            // height of fence
w = 0.1;            // width of fence
t = 5.0;            //legth of tail

fx1 = 0.00;         // front of fence
fx2 = fx1 + w;      // width of fence
fx3 = fx2 + t;      // legth of tail

fy1 = -w/2;         // base of tail
fy2 =  w/2;         // top of tail
fy3 =  fy1-h;       // base of fence
fy4 =  fy2+h;       // base of fence

//------ Domain-------------------------
Point(1) = { xmin,  ymin, zall, lc*1.0};
Point(2) = { xmax,  ymin, zall, lc*1.0};
Point(3) = { xmin,  ymax, zall, lc*1.0};
Point(4) = { xmax,  ymax, zall, lc*1.0};

//---------------------------------------
Point(5) = { fx1,  fy3, zall, 0.025};  // base front
Point(6) = { fx1,  fy4, zall, 0.025};  // top  front

Point(7) = { fx2,  fy4, zall, 0.025};  // top rear 
Point(8) = { fx2,  fy2, zall, 0.050};  // top rear base

Point(9)  = { fx3,  fy2, zall, 0.050};  // tail top 
Point(10) = { fx3,  fy1, zall, 0.050};  // tail base

Point(11) = { fx2,  fy1, zall, 0.025};  // bottom rear base
Point(12) = { fx2,  fy3, zall, 0.050};  // bottom rear

//--------PML----------------------------
Point(13) = { xmin-pmlx1,  ymin, zall, lc*1.0};
Point(14) = { xmin-pmlx1,  ymax, zall, lc*1.0};
Point(15) = { xmax+pmlx2,  ymin, zall, lc*1.0};
Point(16) = { xmax+pmlx2,  ymax, zall, lc*1.0};

Point(17) = { xmin,  ymin-pmly1, zall, lc*1.0};
Point(18) = { xmax,  ymin-pmly1, zall, lc*1.0};
Point(19) = { xmin,  ymax+pmly2, zall, lc*1.0};
Point(20) = { xmax,  ymax+pmly2, zall, lc*1.0};

Point(21) = { xmin-pmlx1,  ymin-pmly1, zall, lc*1.0};
Point(22) = { xmax+pmlx2,  ymin-pmly1, zall, lc*1.0};
Point(23) = { xmin-pmlx1,  ymax+pmly2, zall, lc*1.0};
Point(24) = { xmax+pmlx2,  ymax+pmly2, zall, lc*1.0};

//---------------------------------------
//Domain
Line(1) = {1,2};
Line(2) = {2,4};
Line(3) = {4,3};
Line(4) = {3,1};

//Fence
Line(5)  = { 5, 6};
Line(6)  = { 6, 7};
Line(7)  = { 7, 8};
Line(8)  = { 8, 9};
Line(9)  = { 9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12, 5};

//PML x
Line(13)  = {13,14};
Line(14)  = {15,16};
Line(15)  = { 1,13};
Line(16)  = { 3,14};
Line(17)  = { 2,15};
Line(18)  = { 4,16};

//PML y
Line(19)  = {17,18};
Line(20)  = {19,20};
Line(21)  = { 1,17};
Line(22)  = { 2,18};
Line(23)  = { 3,19};
Line(24)  = { 4,20};

//PML xy
Line(25)  = {13,21};
Line(26)  = {21,17};
Line(27)  = {18,22};
Line(28)  = {22,15};
Line(29)  = {14,23};
Line(30)  = {23,19};
Line(31)  = {20,24};
Line(32)  = {24,16};


//---------------------------------------
// main fluid domain
//---------------------------------------
LLDomain=newll; Line Loop(LLDomain) = {1,2,3,4};
LLFence=newll;  Line Loop(LLFence) = {5,6,7,8,9,10,11,12};
PSDomain=news;  Plane Surface(PSDomain) = {LLDomain, LLFence};

//---------------------------------------
// 12 PML regions
//---------------------------------------
LLpmlx1=newll;      Line Loop(LLpmlx1) = {15, 13, -16, 4};
PSpmlx1=news;   Plane Surface(PSpmlx1) = {LLpmlx1};
LLpmlx2=newll;      Line Loop(LLpmlx2) = {2, 18, -14, -17};
PSpmlx2=news;   Plane Surface(PSpmlx2) = {LLpmlx2};

LLpmly1=newll;      Line Loop(LLpmly1) = {1, 22, -19, -21};
PSpmly1=news;   Plane Surface(PSpmly1) = {LLpmly1};
LLpmly2=newll;      Line Loop(LLpmly2) = {3, 23, 20, -24};
PSpmly2=news;   Plane Surface(PSpmly2) = {LLpmly2};

LLpmlxy1=newll;      Line Loop(LLpmlxy1) = {15, 25, 26, -21};
PSpmlxy1=news;   Plane Surface(PSpmlxy1) = {LLpmlxy1};
LLpmlxy2=newll;      Line Loop(LLpmlxy2) = {22, 27, 28, -17};
PSpmlxy2=news;   Plane Surface(PSpmlxy2) = {LLpmlxy2};
LLpmlxy3=newll;      Line Loop(LLpmlxy3) = {16, 29, 30, -23};
PSpmlxy3=news;   Plane Surface(PSpmlxy3) = {LLpmlxy3};
LLpmlxy4=newll;      Line Loop(LLpmlxy4) = {18, -32, -31, -24};
PSpmlxy4=news;   Plane Surface(PSpmlxy4) = {LLpmlxy4};

Physical Surface("Domian",9) = {PSDomain};
Physical Surface("PMLX",100) = {PSpmlx1,PSpmlx2};
Physical Surface("PMLY",200) = {PSpmly1,PSpmly2};
Physical Surface("PMLXY",300) = {PSpmlxy1,PSpmlxy2,PSpmlxy3,PSpmlxy4};
Physical    Line("Wall", 1) = {5,6,7,8,9,10,11,12};
Physical    Line("Inflow",2) = {29, 13, 25};
Physical    Line("Outflow",3) = {28, 14, 32};
Physical    Line("Slip",4) = {30, 20, 31, 27, 19, 26};

//---------------------------------------------------------
// use Fields to tune wake density
//---------------------------------------------------------

Field[1] = Box;
Field[1].VIn  =  0.035;
Field[1].VOut =  0.50;
Field[1].XMin = -0.8;
Field[1].XMax =  2.5;
Field[1].YMin =  -1.8;
Field[1].YMax =   1.8;

Field[2] = Box;
Field[2].VIn  =  0.05;
Field[2].VOut =  0.50;
Field[2].XMin = -2.2;
Field[2].XMax =  5.5;
Field[2].YMin =  -2.7;
Field[2].YMax =   2.7;

Field[3] = Box;
Field[3].VIn  =  0.05;
Field[3].VOut =  0.50;
Field[3].XMin = -3.3;
Field[3].XMax =  0.0;
Field[3].YMin =  -1.0;
Field[3].YMax =   1.0;

//Field[4] = Box;
//Field[4].VIn  =  0.10;
//Field[4].VOut =  0.50;
//Field[4].XMin = -4.5;
//Field[4].XMax =  8.0;
//Field[4].YMin =  0.0;
//Field[4].YMax =  3.0;

Field[5] = Box;
Field[5].VIn  =  0.15;
Field[5].VOut =  0.50;
Field[5].XMin = -4.0;
Field[5].XMax = 12.0;
Field[5].YMin = -4.0;
Field[5].YMax =  4.0;

Field[6] = Box;
Field[6].VIn  =  0.02;
Field[6].VOut =  0.50;
Field[6].XMin = -0.2;
Field[6].XMax =  0.5;
Field[6].YMin =  0.8;
Field[6].YMax =  1.2;

Field[7] = Box;
Field[7].VIn  =  0.02;
Field[7].VOut =  0.50;
Field[7].XMin = -0.2;
Field[7].XMax =  0.5;
Field[7].YMin =  -1.2;
Field[7].YMax =  -0.8;


Field[8] = Min;
Field[8].FieldsList = {1,2,3,4,5,6,7};
Background Field = 8;


