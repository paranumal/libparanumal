cl__1 = 0.25;

xn  = DefineNumber[4]; 
yn  = DefineNumber[6]; 
Point(1) = {-0.5,-0.5, 0, 1.0};
Point(2) = { 1.0,-0.5, 0, 1.0};
Point(3) = { 1.0, 1.5, 0, 1.0};
Point(4) = {-0.5, 1.5, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line {1,3} = (xn+1) Using Progression 1;
Transfinite Line {2,4} = (yn+1) Using Progression 1;

Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
Transfinite Surface {6};
Recombine Surface {6};

Physical Line("Inflow",2) = {1, 3, 4};
Physical Line("Outflow",3) = {2};
Physical Surface("Domain") = {6};
