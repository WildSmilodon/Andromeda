//
//  gmsh mreza za valj
//

// radij valja
r1=0.5;
r2=1.0;

// dolzina valja
l1=1;
l2=1;

// razdalja med valji
d1=-0.5;
d2=-1.0;

// Velikost elementa
lcar=0.1;


Point(1) = {0,0,d1,lcar};
Point(2) = {r1,0,d1,lcar};
Point(3) = {0,r1,d1,lcar};
Point(4) = {-r1,0,d1,lcar};
Point(5) = {0,-r1,d1,lcar};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

//+
Line Loop(1) = {4, 1, 2, 3};
//+
Surface(1) = {-1};
//+
Rotate {{0, 1, 0}, {0, 0, l1}, Pi} {
  Duplicata { Point{3}; Point{2}; Point{1}; Point{4}; Point{5}; Line{1}; Line{2}; Line{4}; Line{3}; Surface{1}; }
}
//+
Line(9) = {6, 3};
//+
Line(10) = {7, 4};
//+
Line(11) = {2, 9};
//+
Line(12) = {10, 5};
//+
Line Loop(2) = {6, -11, 1, -9};
//+
Surface(10) = {2};
//+
Line Loop(3) = {2, -10, 5, 9};
//+
Surface(11) = {3};
//+
Line Loop(4) = {3, -12, 7, 10};
//+
Surface(12) = {4};
//+
Line Loop(5) = {12, 4, 11, 8};
//+
Surface(13) = {5};
//+
Physical Surface("noter") = {1, 10, 11, 12, 13, 9};


Point(81) = {0,0,d2,lcar};
Point(82) = {r2,0,d2,lcar};
Point(83) = {0,r2,d2,lcar};
Point(84) = {-r2,0,d2,lcar};
Point(85) = {0,-r2,d2,lcar};

Circle(81) = {82,81,83};
Circle(82) = {83,81,84};
Circle(83) = {84,81,85};
Circle(84) = {85,81,82};

//+
Line Loop(81) = {84, 81, 82, 83};
//+
Surface(81) = {-81};
//+
Rotate {{0, 1, 0}, {0, 0, l2}, Pi} {
  Duplicata { Point{83}; Point{82}; Point{81}; Point{84}; Point{85}; Line{81}; Line{82}; Line{84}; Line{83}; Surface{81}; }
}
//+
Line(89) = {86, 83};
//+
Line(90) = {87, 84};
//+
Line(91) = {82, 89};
//+
Line(92) = {90, 85};
//+
Line Loop(82) = {86, -91, 81, -89};
//+
Surface(90) = {82};
//+
Line Loop(83) = {82, -90, 85, 89};
//+
Surface(91) = {83};
//+
Line Loop(84) = {83, -92, 87, 90};
//+
Surface(92) = {84};
//+
Line Loop(85) = {92, 84, 91, 88};
//+
Surface(93) = {85};
//+
Physical Surface("spredaj") = {81};
//+
Physical Surface("stena") = {90, 91, 92, 93};
//+
Physical Surface("zadaj") = {89};
