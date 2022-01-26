# Andromeda, a BEM code

Solves Laplace diffusion and Stokes flow problems with the Boundary Element Method

written by Jure Ravnik, http://kepoi.fs.um.si/en/p/jure-ravnik 


### Papers written using this method:

* A. Šušnjara, O. Verhnjak, D. Poljak, M. Cvetković, J. Ravnik (2021). Stochastic-deterministic boundary element modelling of transcranial electric stimulation using a three layer head model. Engineering Analysis with Boundary Elements, Volume 123, 1 February 2021, Pages 70-83 
doi:10.1016/j.enganabound.2020.11.010

* Anna Šušnjara, Ožbej Verhnjak, Dragan Poljak, Mario Cvetković, Jure Ravnik (2022). Uncertainty quantification and sensitivity analysis of transcranial electric stimulation for 9-subdomain human head model. Engineering Analysis with Boundary Elements (2022), Vol 135, Pages 1-11, doi:10.1016/j.enganabound.2021.10.026

* Jure Ravnik, Mitja Štrakl, Jana Wedel, Paul Sterinamnn, Matjaž Hriberšek: Stokes Flow Induced Drag and Torque on Asbestos Particles can not be estimated by a Simplistic Ellipsoidal Approximation, to appear in 2022

* Mitja Štrakl, Matjaž Hriberšek, Jana Wedel, Paul Sterinamnn, Jure Ravnik: A model for translation and rotation resistance tensors for superellipsoidal particles in Stokes flow, to appear in 2022/2023



### Installation & test:

```
sudo apt-get update
sudo apt install gfortran make libblas-dev liblapack-dev
git clone https://github.com/WildSmilodon/Andromeda.git
cd Andromeda/src
make
cd ../run
./runAll
```


### Paraview is used for results postprocessing

* Install Paraview

```
sudo apt-get update
sudo apt-get install paraview
```

* Open *.vtu files with Paraview to visualise results


### Input files

* ```and.inp``` - main input file
* ```problemName.bic``` - boundary and initial conditions file
* ```problemName.msh``` - mesh file 

### Mesh files

Andromeda uses meshes produced by gmsh https://gmsh.info/

* Install gmsh

```
sudo apt-get update
sudo apt-get install gmsh
```

* Look at examples to study *.geo files, which are an input for gmsh or examine gmsh user manual
* To produce *.msh needed by Andromeda run ```gmsh -2 meshName.geo```

### Getting help

* running Andromeda with ```h``` command line argument displays help.