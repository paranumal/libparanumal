 res = DefineNumber[0.025];
 
 xn = DefineNumber[5];
 yn = DefineNumber[5];
  
 Point(1) = {-1.0, -1.0, 0, res};
 Point(2) = {-1.0, 1.0, 0, res};
 Point(3) = {1.0, 1.0, 0, res};
 Point(4) = {1.0, -1.0, 0, res};
 Point(5) = {-1.0, 0, 0, res};
 Point(6) = {0, -1.0, 0, res};
 Point(7) = {0, 1.0, 0, res};
 Point(8) = {1.0, 0.0, 0, res};
 
 // Point(1) = {-0.5, -0.5, 0, res};
 // Point(2) = {-0.5, 0.5, 0, res};
 // Point(3) = {0.5, 0.5, 0, res};
 // Point(4) = {0.5, -0.5, 0, res};
 // Point(5) = {-0.5, 0, 0, res};
 // Point(6) = {0, -0.5, 0, res};
 // Point(7) = {0, 0.5, 0, res};
 // Point(8) = {0.5, 0.0, 0, res};



 Line(1) = {1, 6};
 Line(2) = {6, 4};
 Line(3) = {4, 8};
 Line(4) = {8, 3};
 Line(5) = {3, 7};
 Line(6) = {7, 2};
 Line(7) = {2, 5};
 Line(8) = {5, 1};
 
 //+
 Transfinite Line {1, 6} = (xn/2+1) Using Progression 1;
 Transfinite Line {2, 5} = (xn/2+1) Using Progression 1;
 Transfinite Line {3, 8} = (yn/2+1) Using Progression 1;
 Transfinite Line {4, 7} = (yn/2+1) Using Progression 1;
 
 Line Loop(9) = {1, 2, 3, 4, 5, 6, 7, 8};
 Plane Surface(9) = {9};
 Transfinite Surface {9} = {1,2,3,4};
 
 Physical Surface("Domain",9) = {9};
 Physical Line("Inflow",2) = {2, 4, 6, 8};
 Physical Line("Outflow",3) = {1, 3, 5, 7};
 // Physical Line("Inflow",2) = {3, 4, 5, 6};
 // Physical Line("Outflow",3) = {1, 2, 7, 8};
 Recombine Surface {9};

