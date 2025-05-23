# Andromeda, a BEM code

Solves Laplace diffusion and Stokes flow problems with the Boundary Element Method

written by Jure Ravnik, http://kepoi.fs.um.si/en/p/jure-ravnik 

### Cite as

Ravnik, J. (2025). Andromeda, a BEM based Stokes flow solver (v1.7). Zenodo. https://doi.org/10.5281/zenodo.14801782

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14801782.svg)](https://doi.org/10.5281/zenodo.14801782)

### Papers written using this method:

* A. Šušnjara, O. Verhnjak, D. Poljak, M. Cvetković, J. Ravnik (2021). Stochastic-deterministic boundary element modelling of transcranial electric stimulation using a three layer head model. Engineering Analysis with Boundary Elements, Volume 123, 1 February 2021, Pages 70-83 
doi:10.1016/j.enganabound.2020.11.010

* Anna Šušnjara, Ožbej Verhnjak, Dragan Poljak, Mario Cvetković, Jure Ravnik (2022). Uncertainty quantification and sensitivity analysis of transcranial electric stimulation for 9-subdomain human head model. Engineering Analysis with Boundary Elements (2022), Vol 135, Pages 1-11, doi:10.1016/j.enganabound.2021.10.026

* J. Ravnik, M. Štrakl, J. Wedel, M. Steinmann, M. Hriberšek. (2022). Stokes flow induced drag and torque on asbestos like fibres cannot be estimated by a simplistic ellipsoidal approximation WIT Transactions on Engineering Sciences, vol. 134, pp. 33 - 44; doi:10.2495/BE450031

* Štrakl, M., Hriberšek, M., Wedel, J., Steinmann, P., Ravnik, J. (2022). A Model for Translation and Rotation Resistance Tensors for Superellipsoidal Particles in Stokes Flow. J. Mar. Sci. Eng. 2022, 10, 369. doi:10.3390/jmse10030369

* J. Ravnik (2023). Analytical expressions for singular integrals arising from the 3D Laplace and Stokes kernels when using constant or linear triangular and quadrilateral boundary elements. Engineering Analysis with Boundary Elements, Volume 154, September 2023, Pages 47-53. doi:10.1016/j.enganabound.2023.02.057

### Uses LIS, blas and mpich libraries, compiles under gfortran

```
sudo apt-get update
sudo apt install gfortran make libblas-dev liblapack-dev mpich
```

Download and install LIS library from: https://www.ssisc.org/lis

LIS library is used with modified header file, see ```./inc/lisf.h```


### Installation & test:

```
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

* running Andromeda with ```-h``` command line argument displays help.

### Docker

To make a docker image of Andromeda, using the Dockerfile provided, run

Build

```docker build -t andromeda .```

Run

```docker run -it andromeda```

To have access to a folder from inside docker

```docker run -v $(pwd)/my_files:/workspace/data -it andromeda```
