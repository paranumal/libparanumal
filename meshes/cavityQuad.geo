res = DefineNumber[0.025];
 
 xn = DefineNumber[32];
 yn = DefineNumber[32];
  
 Point(1) = {-1.0, -1.0, 0, res};
 Point(2) = {-1.0, 1.0, 0, res};
 Point(3) = {1.0, 1.0, 0, res};
 Point(4) = {1.0, -1.0, 0, res};
  

 Line(1) = {1, 2};
 Line(2) = {2, 3};
 Line(3) = {3, 4};
 Line(4) = {4, 1};
 
 //+
 Transfinite Line {1, 3} = (xn+1) Using Progression 1;
 Transfinite Line {2, 4} = (yn+1) Using Progression 1;
 
 Line Loop(9) = {1, 2, 3, 4};
 Plane Surface(9) = {9};
 Transfinite Surface {9} = {1,2,3,4};
 
 Physical Surface("Domain",9) = {9};
 // Physical Line("Inflow",1) = {1,2,3,4};
 Physical Line("Inflow",1) = {1,4};
 Physical Line("Outflow",2) = {2,3};
 Recombine Surface {9};

