//
//  gmsh mreza za kroglo v valju
//

// Novi gmsh ne zasuka normal, uporabi ali
// gmsh -2 -format msh2 sphere-in-a-cyl.geo
// ali
// ukaz ReverseMesh Surface{16}; to naredi

// radij valja
r2=1.5;

// Dolzina valja
d = 6;

// Velikost elementa v valju
lcar=0.25;


s=0.05; // Velikost elementa v krogli
Mesh.CharacteristicLengthFromCurvature = 0.05;
r=0.5; // radij krogle



// *********************************************************************
// **************** KROGLA *********************************************
// *********************************************************************

Point(1) = {0, 0, 0, s};
Point(2) = {-r, 0, 0, s};
Point(3) = {r, 0, 0, s};
Point(4) = {0, -r, 0, s};
Point(5) = {0, r, 0, s};
Point(6) = {0, 0, r, s};
Point(7) = {0, 0, -r, s};
Circle(1) = {3, 1, 5};
Circle(2) = {2, 1, 5};
Circle(3) = {2, 1, 4};
Circle(4) = {4, 1, 3};
Circle(5) = {6, 1, 3};
Circle(6) = {6, 1, 2};
Circle(7) = {2, 1, 7};
Circle(8) = {7, 1, 3};
Circle(9) = {5, 1, 6};
Circle(10) = {6, 1, 4};
Circle(11) = {4, 1, 7};
Circle(12) = {7, 1, 5};
Line Loop(13) = {9, 5, 1};
Surface(14) = {13};
Line Loop(15) = {2, 9, 6};
Surface(16) = {15};
Line Loop(17) = {7, 12, -2};
Surface(18) = {17};
Line Loop(19) = {12, -1, -8};
Surface(20) = {19};
Line Loop(21) = {6, 3, -10};
Surface(22) = {21};
Line Loop(23) = {10, 4, -5};
Surface(24) = {23};
Line Loop(25) = {11, -7, 3};
Surface(26) = {25};
Line Loop(27) = {8, -4, 11};
Surface(28) = {27};

// Reverse normals
ReverseMesh Surface{16}; // this was done using "-" in Physical Sruface in v2 mesh, in v4 it is not
ReverseMesh Surface{18}; // this was automatic in v2 mesh, in v4 it is not
ReverseMesh Surface{26}; // this was automatic in v2 mesh, in v4 it is not

//Physical surfaces of domian
Physical Surface("krogla") = {14, 16, 18, 20, 22, 24, 26, 28};


// *********************************************************************
// **************** VALJ ***********************************************
// *********************************************************************


// začetek valja na z osi
d2=-d/2;
// tocka okoli katere zasukamo osnovno ploskev (ki je pri z=d2)
l2=0;


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
