//
//
// number of layers
n=5;
// length
l=1;

Point(1) = {0,0,0,0.2};
Extrude {l,0,0} {
  Point{1}; Layers{n};
}
Extrude {0,l,0} {
  Line{1}; Layers{n}; Recombine;
}
Extrude {0,0,l} {
  Surface{5}; Layers{n}; Recombine;
} 
Physical Surface("spodaj") = {-5};
Physical Surface("spredaj") = {14};
Physical Surface("desno") = {18};
Physical Surface("zadaj") = {22};
Physical Surface("levo") = {26};
Physical Surface("zgoraj") = {27};


