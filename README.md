<p align="center">
  <img src="https://github.com/user-attachments/assets/4f06bbdf-ef2f-4500-b775-fd8cfe9dd590" width=250 />

</p>
<h3 align="center">Growth Balance Analysis for C++</h3>

<p align="center">
<br />
<a href="https://github.com/charlesrocabert/gbacpp/releases/latest"><img src="https://img.shields.io/github/release/charlesrocabert/gbacpp/all.svg" /></a>&nbsp;
<a href="https://github.com/charlesrocabert/gbacpp/actions/workflows/cmake-multi-platform.yml"><img src="https://github.com/charlesrocabert/gbacpp/actions/workflows/cmake-multi-platform.yml/badge.svg" /></a>&nbsp;
<a href="https://github.com/charlesrocabert/gbacpp/LICENSE.html"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" /></a>
</p>

<p align="center">
  <a href="https://www.cs.hhu.de/en/research-groups/computational-cell-biology" target="_blank"><img src="https://github.com/user-attachments/assets/4e4b3b79-0d6a-4328-9c3f-3497401887e4" width=150 /></a>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <a href="https://www.hhu.de/en/" target="_blank"><img src="https://github.com/user-attachments/assets/7db5c8f7-e37a-415f-88c3-1b06a49e1f28" width=150 /></a>
</p>

-----------------

<p align="justify">
<strong>gbacpp</strong> is a C++ implementation of the growth balance analysis mathematical formalism (GBA; <a href="https://doi.org/10.1371/journal.pcbi.1011156" target="_blank">Dourado et al. 2023</a>).
The software has been optimized to solve large-scale cell growth models (CGMs), for which available solvers for non-linear constraint-based problems usually struggle.
The optimization process relies on a gradient ascent approach, and is preferred for models offering a convex solution space (typically, when all reactions in the metabolic network are linearly independent).
</p>

<p align="justify">
⚠️ Note that CGMs must comply to a standardized format. A tutorial is available in the section <a href="#cgm_format_tutorial">CGM format tutorial</a>.
</p>

<p align="justify">
<strong>gbacpp</strong> will integrate CGM evolutionary algorithms in the near future (see <a href="#roadmap">Roadmap</a>).
</p>

# Table of contents
- [1) Roadmap](#roadmap)
- [2) Installation instructions](#installation_instructions)
  - [2.1) Supported platforms](#supported_platforms)
  - [2.2) Dependencies](#dependencies)
  - [2.3) Software compilation](#software_compilation)
- [3) First usage](#first_usage)
  - [3.1) Find an optimum](#find_optimum)
  - [3.2) Optimization parameters](#optimization_parameters)
  - [3.3) Usage example](#usage_example)
  - [3.4) Ready-to-use examples](#examples)
- [4) CGM format tutorial](#cgm_format_tutorial)
- [5) Copyright](#copyright)
- [6) License](#license)

# 1) Roadmap <a name="roadmap"></a>

| Task | Status |
|---|---|
| Gradient ascent (best for full column-rank CGMs) | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-Done-green.svg"/></a> |
| MCMC algorithm | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-To do-red.svg"/></a> |
| Forward-in-time population level simulations | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-To do-red.svg"/></a> |
| Lineage tracking | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-To do-red.svg"/></a> |

# 2) Installation instructions <a name="installation_instructions"></a>
Download the <a href="https://github.com/charlesrocabert/gbacpp/releases/latest">latest release</a> of <strong>gbacpp</strong>, and save it to a directory of your choice. Open a terminal and use the <code>cd</code> command to navigate to this directory. Then follow the steps below to compile and build the executables.

### 2.1) Supported platforms <a name="supported_platforms"></a>
<strong>gbacpp</strong> software has been primilary developed for Unix/Linux and macOS systems.

### 2.2) Dependencies <a name="dependencies"></a>
* A C++ compiler (GCC, LLVM, ...; C++17 required),
* CMake $\geq$ 3.5 (command line version),
* GSL $\geq$ 2.8 (https://www.gnu.org/software/gsl/).

### 2.3) Software compilation <a name="software_compilation"></a>

#### • User mode
To compile <strong>gbacpp</strong>, run the following instructions on the command line:

    cd cmake/

and

    bash make.sh

#### • Debug mode
To compile the software in DEBUG mode, use <code>make_debug.sh</code> script instead of <code>make.sh</code>:

    bash make_debug.sh

This mode should only be used for test or development phases.

#### • Executable files emplacement
Binary executable files are in <code>build/bin</code> folder.

#### • Delete compiled files
To clean compiled files and binary executables, run:

    bash make_clean.sh

# 3) First usage <a name="first_usage"></a>
Once <strong>gbacpp</strong> has been compiled, follow the next steps for a first usage of the software.

### 3.1) Find an optimum <a name="find_optimum"></a>
To run a gradient ascent optimization on a CGM, execute the following command line:

    ./build/bin/find_optimum <parameters>

The command line parameters are described below. The description is also available by executing the following command line in a terminal:

    ./build/bin/find_optimum -h

### 3.2) Optimization parameters <a name="optimization_parameters"></a>

- <code>-h</code>, <code>--help</code>: Print the help, then exit,
- <code>-v</code>, <code>--version</code>: Print the current version, then exit,

#### • Mandatory parameters
- <code>-path</code>, <code>--model-path</code>: Specify the path of the CGM to be loaded,
- <code>-name</code>, <code>--model-name</code>: Specify the name of the CGM to be loaded,
- <code>-condition</code>, <code>--condition</code>: Specify the external condition. If <code>all</code> is selected, all conditions are applied.

#### • Optional parameters
- <code>-print</code>, <code>--print-optimum</code>: Indicates if the optimum should be printed in the standard output. This option is useful to pass the result to another program (<code>-verbose</code> option should not be used),
- <code>-write</code>, <code>--write-trajectory</code>: Indicates if the trajectory should be written as output files. Tracking the optimization trajectory can be useful during tests,
- <code>-output</code>, <code>--output-path</code>: Specify the path of output files,
- <code>-tol</code>, <code>--tolerance</code>: Specify the tolerance value ($10^{-10}$ by default),
- <code>-stable-count</code>, <code>--stable-count</code>: Specify the maximal number of iterations with unchanged mu ($10,000$ by default),
- <code>-maxt</code>, <code>--max-time</code>: Specify the maximal trajectory time ($100,000$ by default),
- <code>-verbose</code>, <code>--verbose</code>: Indicates if the program should run in verbose mode (can conflict with the option <code>-print</code>).

### 3.3) Usage example <a name="usage_example"></a>

In this example, growth rate optimums are calculated for all the external conditions of the <em>Escherichia coli</em> toy model EC12b (see folder <code>./examples/toy_models</code>).

First, navigate to the folder <code>./examples</code> using the <code>cd</code> command:

    cd ./examples

Then, call the optimization algorithm:

    ../build/bin/find_optimum -path ./toy_models -name EC12b -condition all -output ./output -verbose

Here, optimums are calculated for all conditions, and saved in the folder <code>./examples/output</code>. Verbose mode is activated to get insights in the optimization process.

### 3.4) Ready-to-use examples <a name="examples"></a>
Ready-to-use examples are available in the folder <code>./examples</code> (navigate to the folder <code>examples</code> using the <code>cd</code> command):

• <code>model_A_condition_1.sh</code>: This script will run a single gradient ascent on model A in external condition 1 (2 reactions, 2 metabolites). You can execute it using the following command line:

    bash model_A_condition_1.sh

At the end of the optimization, CSV files are written in the folder <code>./examples/output</code>. You can edit the parameter values at will to test the behaviour of the gradient ascent. See below for a full description of the parameters.

• <code>model_B_all_conditions.sh</code>: This script will calculate the optimal growth rate on model B for all external conditions. You can execute it using the following command line:

    bash model_B_all_conditions.sh

All the optimums are written in the folder <code>./examples/output</code>.

## 4) CGM format tutorial <a name="cgm_format_tutorial"></a>

Coming soon. You can refer to the https://cellgrowthsim.com/ tutorial.

## 5) Copyright <a name="copyright"></a>
Copyright © 2024-2025 Charles Rocabert. All rights reserved.

## 6) License <a name="license"></a>
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
