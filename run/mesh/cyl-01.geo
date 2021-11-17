//
//  gmsh mreza za valj
//

// radij valja
r=0.5;

// dolzina valja
l=1;

// Velikost elementa
lcar=0.1;

Point(1) = {0,0,0,lcar};
Point(2) = {r,0,0,lcar};
Point(3) = {0,r,0,lcar};
Point(4) = {-r,0,0,lcar};
Point(5) = {0,-r,0,lcar};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Extrude {0,0,l} {
  Surface{6};
}
Physical Surface("spredaj") = {-6};  // minus zato, da so normale prav
Physical Surface("stena") = {27, 15, 23, 19};
Physical Surface("zadaj") = {28};

// Izracuna integral po povrsini
//Mesh 2; // mesh
//Plugin(NewView).Run; // create new post-pro view
//Plugin(ModifyComponent).View = 0;
//Plugin(ModifyComponent).Expression = "1"; 
//Plugin(ModifyComponent).Run; // assign 1 on this view
//Plugin(Integrate).Dimension = 2;
//Plugin(Integrate).Run; // integrate 1 (only on 2D elements)
//Printf("Area = %g", View[1].Max); // print out the result
