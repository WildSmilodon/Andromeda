# Andromeda, a BEM code


written by Jure Ravnik, http://kepoi.fs.um.si/en/p/jure-ravnik 


Papers written using this method:

* A. Šušnjara, O. Verhnjak, D. Poljak, M. Cvetković, J. Ravnik (2021). Stochastic-deterministic boundary element modelling of transcranial electric stimulation using a three layer head model. Engineering Analysis with Boundary Elements, Volume 123, 1 February 2021, Pages 70-83 
doi:10.1016/j.enganabound.2020.11.010

* Anna Šušnjara, Ožbej Verhnjak, Dragan Poljak, Mario Cvetković, Jure Ravnik (2022). Uncertainty quantification and sensitivity analysis of transcranial electric stimulation for 9-subdomain human head model. Engineering Analysis with Boundary Elements (2022), Vol 135, Pages 1-11, doi:10.1016/j.enganabound.2021.10.026

Installation & test:

sudo apt install gfortran make libblas-dev liblapack-dev<br>
git clone https://github.com/WildSmilodon/Andromeda.git<br>
cd Andromeda/src<br>
make<br>
cd ../run<br>
./runAll<br>
