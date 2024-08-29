<h1 align="center">&sigma;FGM</h1>
<p align="center">
<img src="logo/logo.png" width="150">
<br/>
<em>Individual-based implementation of Fisher's geometric model with evolvable phenotypic noise</em>
<br/><br/>
<a href="https://github.com/charlesrocabert/SigmaFGM/releases/latest"><img src="https://img.shields.io/github/release/charlesrocabert/SigmaFGM/all.svg" /></a>&nbsp;
<a href="https://github.com/charlesrocabert/SigmaFGM/actions"><img src="https://github.com/charlesrocabert/SigmaFGM/workflows/CMake/badge.svg" /></a>&nbsp;
<a href="https://github.com/charlesrocabert/SigmaFGM/LICENSE.html"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" /></a>
</p>

-----------------

<p align="justify">
Experimental studies demonstrate the existence of phenotypic diversity despite constant genotype and environment. Theoretical models based on a single phenotypic character predict that during an adaptation event, phenotypic noise should be positively selected far from the fitness optimum because it increases the fitness of the genotype, and then be selected against when the population reaches the optimum. It is suggested that because of this fitness gain, phenotypic noise should promote adaptive evolution. However, it is unclear how the selective advantage of phenotypic noise is linked to the rate of evolution, and whether any advantage would hold for more realistic, multi-dimensional phenotypes. Indeed, complex organisms suffer a cost of complexity, where beneficial mutations become rarer as the number of phenotypic characters increases (<a href="https://doi.org/10.1111/evo.14083">Rocabert et al. 2020</a>).
</p>

<p align="justify">
&sigma;FGM simulates adaptive evolution in Fisher's geometric model with an evolvable phenotypic noise. A population of individuals is placed under stabilizing selection and must evolve towards a fitness optimum. The fitness function is configurable and can adopt non-Gaussian shapes. Phenotypic noise is modeled by an evolvable multivariate normal distribution. Simulations are fully configurable. See below for <a href="#installation">installation instructions</a> and a <a href="#first_usage">first usage</a>, and <a href="https://doi.org/10.1111/evo.14083">Rocabert et al. (2020)</a> for a full description of the underlying mathematical model.
</p>

<p align="justify">
&sigma;FGM has been developed by Charles Rocabert, Guillaume Beslon, Carole Knibbe and Samuel Bernard.
</p>

## Table of contents
- [Publications](#publications)
- [Installation instructions](#installation)
- [First usage](#first_usage)
- [Copyright](#copyright)
- [License](#license)

## Publications <a name="publications"></a>

• Rocabert C., Beslon G., Knibbe C. and Bernard S. (2020). Phenotypic Noise and the Cost of Complexity, _Evolution_, in press. https://doi.org/10.1111/evo.14083

## Installation instructions <a name="installation"></a>

Download the <a href="https://github.com/charlesrocabert/SigmaFGM/releases/latest">latest release</a> of &sigma;FGM, and save it to a directory of your choice. Open a terminal and use the <code>cd</code> command to navigate to this directory. Then follow the steps below to compile and build the executables.

### 1. Supported platforms
&sigma;FGM software has been developed for Unix/Linux and macOS systems.

### 2. Required dependencies
* A C++ compiler (GCC, LLVM, ...)
* CMake (command line version)
* GSL
* CBLAS

### 3. Software compilation
Download the latest release of &sigma;FGM, and save it to a directory of your choice. Open a terminal and use the <code>cd</code> command to navigate to this directory. Then follow the steps below to compile and build the executables.

#### User mode
To compile &sigma;FGM, run the following instructions on the command line:

    cd cmake/

and

    bash make.sh

#### Debug mode
To compile the software in DEBUG mode, use <code>make_debug.sh</code> script instead of <code>make.sh</code>:

    bash make_debug.sh

This mode should only be used for test or development phases.

#### Executable files emplacement
Binary executable files are in <code>build/bin</code> folder.

## First usage <a name="first_usage"></a>
Open a terminal and use the <code>cd</code> command to navigate to the <code>example</code> directory. Then follow the steps below for a first usage.

As a first step, one can execute the Python file available in the <code>example</code> directory:

    python3 run_SigmaFGM_simulation.py

The parameters of the simulation can be edited in this file (see below for a description of the parameters).

To directly run a simulation from the main executable, use the following command line:

    ../build/bin/SigmaFGM_simulation <parameters>

The command line parameters are described below (see <a href="https://doi.org/10.1111/evo.14083">Rocabert et al. 2020</a> for details on the mathematical terms). Other non-mandatory parameters are also available, and are described by writing the following command line in a terminal:

    ../build/bin/SigmaFGM_simulation -h

#### Parameters:
- <code>-h</code>, <code>--help</code>: Print this help, then exit
- <code>-v</code>, <code>--version</code>: Print the current version, then exit
- <code>-seed</code>, <code>--seed</code>: Specify the PRNG seed (**mandatory**, random if 0)
- <code>-g</code>, <code>--generations</code>: Specify the number of generations (**mandatory**)
- <code>-nbdim</code>, <code>--nb-dimensions</code>: Specify the number of dimensions _n_ (**mandatory**)
- <code>-alpha</code>, <code>--alpha</code>: Specify the &alpha; parameter of the fitness function (**mandatory**)
- <code>-beta</code>, <code>--beta</code>: Specify the &beta; parameter of the fitness function (**mandatory**)
- <code>-Q</code>, <code>--Q</code>: Specify the _Q_ parameter of the fitness function (**mandatory**)
- <code>-popsize</code>, <code>--population-size</code>: Specify the population size _N_ (**mandatory**)
- <code>-initmu</code>, <code>--initial-mu</code>: Specify the initial distance of the mean phenotype **&mu;** from the fitness optimum (**mandatory**)
- <code>-initsigma</code>, <code>--initial-sigma</code>: Specify the initial phenotypic noise amplitude (**mandatory**)
- <code>-inittheta</code>, <code>--initial-theta</code>: Specify the initial phenotypic noise orientation (**mandatory**, usually 0)
- <code>-mmu</code>, <code>--m-mu</code>: Specify **&mu;** mutation rate (**mandatory**)
- <code>-msigma</code>, <code>--m-sigma</code>: Specify **&sigma;** mutation rate (**mandatory**)
- <code>-mtheta</code>, <code>--m-theta</code>: Specify **&theta;** mutation rate (**mandatory**)
- <code>-smu</code>, <code>--s-mu</code>: Specify **&mu;** mutation size (**mandatory**)
- <code>-ssigma</code>, <code>--s-sigma</code>: Specify **&sigma;** mutation size (**mandatory**)
- <code>-stheta</code>, <code>--stheta</code>: Specify **&theta;** mutation size (**mandatory**)
- <code>-noise</code>, <code>--noise-type</code>: Specify the type of phenotypic noise (**mandatory**, NONE/ISOTROPIC/UNCORRELATED/FULL)

Note that setting <code>-noise</code> to <code>NONE</code> leads to a simulation with the classical Fisher's geometric model.

The software outputs two statistics files during the course of the simulation, containing the mean (<code>mean.txt</code>) and the standard deviation (<code>sd.txt</code>) of some metrics allowing to track the state of the evolving population (see <a href="https://doi.org/10.1111/evo.14083">Rocabert et al. 2020</a> for a full description):
- <code>g</code>: Current generation,
- <code>dmu</code>: Distance of the mean phenotype &mu; from the optimum,
- <code>dz</code>: Distance of the phenotype _z_ from the optimum,
- <code>Wmu</code>: Absolute fitness of the mean phenotype &mu;,
- <code>Wz</code>: Absolute fitness of the phenotype _z_,
- <code>EV</code>: Maximal eigenvalue of the phenotypic noise distribution,
- <code>EV_contrib</code>: Contribution of the maximal eigenvalue to the sum of eigenvalues,
- <code>EV_dot_product</code>: Alignment of the eigenvector corresponding to the maximum eigenvalue with the dirction of the fitness optimum, 
- <code>r_mu</code>: Mutation size on the mean phenotype **&mu;**,
- <code>r_sigma</code>: Mutation size on the phenotypic noise amplitudes **&sigma;**,
- <code>r_theta</code>: Mutation size on the phenotypic rotation angles **&theta;**.

## Copyright <a name="copyright"></a>
Copyright &copy; 2016-2020 Charles Rocabert, Guillaume Beslon, Carole Knibbe, Samuel Bernard.
All rights reserved.

## License <a name="license"></a>
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

