<p align="center">
  <a href="http://mathlab.github.io/ITHACA-DG/" target="_blank" >
    <img alt="ITHACA-DG" src="./docs/logo/ithaca-dg_small.png" width="200" />
  </a>
</p>

## ITHACA-DG - In real Time Highly Advanced Computational Applications for Discontinuous Galerkin - ROMs for HopeFOAM ##

<p align="center">
    <a href="https://www.gnu.org/licenses/lgpl-3.0" target="_blank">
        <img alt="Software License" src="https://img.shields.io/badge/License-LGPL%20v3-blue.svg">
    </a>
</p>

### 0. Introduction
**ITHACA-DG** is an implementation in **HopeFOAM** of several reduced order modelling techniques. **ITHACA-DG** is designed for [**HopeFOAM 0.1**](https://github.com/HopeFOAM/HopeFOAM), which is based on OpenFOAM 4.0  [**OpenFOAM 4.0**](https://openfoam.org/version/4) . 


Linear and non-linear algebra operations which are not already implemented in OpenFOAM are performed with the external library [**Eigen**](http://eigen.tuxfamily.org/index.php?title=Main_Page). The source code of Eigen 3.3.4 is provided together with ITHACA-DG and is located in the [src/thirdyparty/Eigen](./src/thirdparty/Eigen) folder.  For the EigenValue decomposition it is also possible to rely on the [**Spectra-0.6.1**](https://spectralib.org/) library and the source code is provided in the [src/thirdyparty/spectra-0.6.1](./src//thirdparty/spectra-0.6.1) folder.

**ITHACA-DG** has been tested on ubuntu 18.04 but can be easily compiled on any linux distribution with a compiled version of OpenFOAM 4.0. 

### 1. Prerequisites
**ITHACA-DG** requires
* [**HopeFOAM 0.1**](https://github.com/HopeFOAM/HopeFOAM) 


### 2. Installation and usage
First of all you need to source the bashrc file of your installation of **HopeFOAM 0.1**. This is of course depending on the location of your OpenFOAM installation and of your particular version of OpenFOAM
```
source $HOME/HopeFOAM/HopeFOAM-0.1/etc/bashrc
``` 
Then navigate to the folder where you want to install ITHACA-DG such as, for example, the utilities folder of OpenFOAM
```
cd ${FOAM_APP}/utilities
``` 
Now you can clone the **ITHACA-DG** repository inside the selected folder
```
git clone https://github.com/mathLab/ITHACA-DG
```
and you can compile **ITHACA-DG** by navigating inside the src folder and compiling using wmake
```
cd ITHACA-DG
./Allwmake 
```
For a brief description of the classes and methods, you can check the official ITHACA-DG doxygen [documentation](https://mathlab.github.io/ITHACA-DG/).


### 3. [Tutorials](https://mathlab.github.io/ITHACA-DG//examples.html)
Tutorials are provided the [**tutorials** subfolder](./tutorials).
* [**Tutorial 1**] In this tutorial is implemented the development of a parametrized POD-Galerkin method for an unsteady Navier-Stokes problem. The parametrization is on the viscosity. The OpenFOAM full order problem is based on **dgChorin** solver.


### 4. Authors and contributors
**ITHACA-DG** is currently developed and mantained at [SISSA mathLab](http://mathlab.sissa.it/) by [Dr. Andrea Lario](mailto:alario@sissa.it) under the supervision of [Prof. Gianluigi Rozza](mailto:gianluigi.rozza@sissa.it)

Contact us by email for further information or questions about **ITHACA-DG**, or open an ''Issue'' on this website. **ITHACA-DG** is at an early development stage, so contributions improving either the code or the documentation are welcome, both as patches or merge requests on this website.

### 5. License
**ITHACA-DG** is freely available under the GNU LGPL, version 3.

![Google Analytics](https://ga-beacon.appspot.com/UA-66224794-1/rbnics/readme?pixel)
