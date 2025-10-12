<p align="center">
  <img src="https://github.com/user-attachments/assets/2a888a9f-1653-469a-9739-55f4cd9e3b21" width=200 />

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
The software has been optimized to solve large-scale self-replicating cell (SRC) models, for which available solvers for non-linear constraint-based problems usually struggle.
The optimization process relies on a gradient ascent approach, and is preferred for models offering a convex solution space.
</p>

<p align="justify">
:ballot_box_with_check: Note that SRC models must comply to a standardized format. Guidelines are available in the <a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/src_format_tutorial.md" target="_blank">SRC model format tutorial</a>.

:ballot_box_with_check: When building a SRC model, stoichiometric coefficients, and kinetic parameters must be converted following GBA formalism. See the <a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/units_conversion_tutorial.ipynb" target="_blank">units conversion tutorial</a>.

:ballot_box_with_check: The gradient ascent algorithm will be detailed in the documentation soon.

:ballot_box_with_check: <strong>gbacpp</strong> will integrate evolutionary algorithms in the near future (see <a href="#roadmap">Roadmap</a>).

:ballot_box_with_check: A complete tutorial on GBA formalism is also available at https://cellgrowthsim.com/.
</p>

# Table of contents
- [1) Roadmap](#roadmap)
- [2) Contributing](#contributing)
- [3) Installation instructions](#installation_instructions)
  - [3.1) Supported platforms](#supported_platforms)
  - [3.2) Dependencies](#dependencies)
  - [3.3) Installation](#installation)
  - [3.4) Manual software compilation](#manual_software_compilation)
- [4) First usage](#first_usage)
  - [4.1) Why using a gradient ascent](#gradient_ascent)
  - [4.2) Code optimization](#optimization)
  - [4.3) Find an optimum](#find_optimum)
  - [4.4) Optimization parameters](#optimization_parameters)
- [5) SRC model format tutorial](#src_model_format_tutorial)
- [6) Units conversion tutorial](#units_conversion_tutorial)
- [7) Copyright](#copyright)
- [8) License](#license)

# 1) Roadmap <a name="roadmap"></a>

| Task | Status |
|---|---|
| Gradient ascent (best for full column-rank SRC models with minimal support) | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-Done-green.svg"/></a> |
| Reload/restart a gradient ascent | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-Done-green.svg"/></a> |
| Handling ODS models | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-To do-red.svg"/></a> |
| MCMC algorithm | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-To do-red.svg"/></a> |
| Forward-in-time population level simulations | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-To do-red.svg"/></a> |
| Lineage tracking | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-To do-red.svg"/></a> |

# 2) Contributing <a name="contributing"></a>

If you wish to contribute, do not hesitate to reach <a href="mailto:charles DOT rocabert AT hhu DOT de">the developer</href>.

# 3) Installation instructions <a name="installation_instructions"></a>
Download the <a href="https://github.com/charlesrocabert/gbacpp/releases/latest">latest release</a> of <strong>gbacpp</strong>, and save it to a directory of your choice. Open a terminal and use the <code>cd</code> command to navigate to this directory. Then follow the steps below to compile and build the executables.

### 3.1) Supported platforms <a name="supported_platforms"></a>
<strong>gbacpp</strong> software has been primilary developed for Unix/Linux and macOS systems.

### 3.2) Dependencies <a name="dependencies"></a>
* A C++ compiler (GCC, LLVM, ...; C++17 required),
* CMake $\geq$ 3.5 (command line version),
* GSL $\geq$ 2.8 (https://www.gnu.org/software/gsl/).

### 3.3) Installation <a name="installation"></a>

<p align="justify">
Download the <a href="https://github.com/charlesrocabert/gbacpp/releases/latest">latest release</a> of <strong>gbacpp</strong>, and save it into a directory of your choice. Open a terminal and use the <code>cd</code> command to navigate to this directory. To install <strong>gbacpp</strong>, simply call the script <code>install.sh</code> on the command line:
</p>

```
sh install.sh
```

<p align="justify">
The script will compile and install the software in the appropriate folder (usually in a <code>bin</code> directory, such as <code>/usr/local/bin</code>).
The binary executable <code>find_model_optimum</code> should be available in the terminal. If not, you may need to export the binary path. For example:
</p>

```
export PATH="/usr/bin:$PATH"
```

<p align="justify">
:warning: If you want to simply compile the software without installing it into your system, follow the next instructions.
</p>

### 3.4) Manual software compilation <a name="manual_software_compilation"></a>

#### • User mode
To manually compile <strong>gbacpp</strong>, run the following instructions on the command line:

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

# 4) First usage <a name="first_usage"></a>
Once <strong>gbacpp</strong> has been compiled and installed, follow the next steps for a first usage of the software.

### 4.1) Why using a gradient ascent <a name="gradient_ascent"></a>

SRC models are non-linear constraint-based problems. Solving such tasks in much more difficult than linear problems (<em>e.g.</em> FBA problems).
A few solvers are available, and can usually handle small SRC models. However, we noticed that these solvers tend to fail when solving larger GBA problems, and are unable to handle genome-scale models.

As GBA formalism provides analytical solutions to calculate the growth rate gradient, we have implemented <strong>gbacpp</strong> to circumvent this limitation, allowing us to solve large, or even genome-scale SRC models in a reasonable timeframe.

<!--
In <strong>gbacpp</strong>, gradient ascent is timestep-adaptive, ensuring solutions stay consistent. In particular, we have implemented controls when concentrations or fluxes converge to zero to avoid algorithmic traps.
-->

### 4.2) Code optimization <a name="optimization"></a>

The gradient ascent implemented here relies on $\dfrac{\partial \mu}{\partial f}$, the growth rate derivative against the flux fraction vector $f$. As presented in <a href="https://doi.org/10.1371/journal.pcbi.1011156" target="_blank">Dourado et al. (2023)</a>, analytical expressions are available to explicitely calculate these values, while it requires heavy linear algebra.

<strong>gbacpp</strong> combines two benefits of using C++. Compilers natively optimize calculations (<em>e.g.</em>, using vectorization), and we could strongly optimize calculations and memory management.
This approach makes <strong>gbacpp</strong> a fast solution.

### 4.3) Find an optimum <a name="find_optimum"></a>
To run a gradient ascent optimization on a SRC model, execute the following command line:

    find_model_optimum <parameters>

The command line parameters are described below. The description is also available by executing the following command line in a terminal:

    find_model_optimum -h

### 4.4) Optimization parameters <a name="optimization_parameters"></a>

- <code>-h</code>, <code>--help</code>: Print the help, then exit,
- <code>-version</code>, <code>--version</code>: Print the current version, then exit,

#### • Mandatory parameters
- <code>-path</code>, <code>--model-path</code>: Specify the path of the SRC model to be loaded,
- <code>-name</code>, <code>--model-name</code>: Specify the name of the SRC model to be loaded,
- <code>-condition</code>, <code>--condition</code>: Specify the external condition. If <code>all</code> is selected, all conditions are applied.

#### • Optional parameters
- <code>-print</code>, <code>--print-optimum</code>: Indicates if the optimum should be printed in the standard output. This option is useful to pass the result to another program (⚠️ <code>-verbose</code> option should not be used),
- <code>-optimum</code>, <code>--write-optimum</code>: Indicates if the optimum should be written as output files,
- <code>-trajectory</code>, <code>--write-trajectory</code>: Indicates if the trajectory should be written as output files. Tracking the optimization trajectory can be useful during tests,
- <code>-output</code>, <code>--output-path</code>: Specify the path of output files,
- <code>-tol</code>, <code>--tolerance</code>: Specify the tolerance value ($10^{-10}$ by default),
- <code>-mutol</code>, <code>--mu-tolerance</code>: Specify the relative growth rate difference tolerance value ($10^{-10}$ by default),
- <code>-stable</code>, <code>--stable-count</code>: Specify the maximal number of iterations with unchanged growth rate as a stop criterium ($10,000$ by default),
- <code>-max</code>, <code>--max-iter</code>: Specify the maximal number of iterations as a stop criterium ($100,000,000$ by default),
- <code>-reload</code>, <code>--reload</code>: Indicates if the last trajectory point should be used as q0
- <code>-restart</code>, <code>--restart</code>: Indicates if the last trajectory point should be used as a fresh start
- <code>-previous</code>, <code>--use-previous-sol</code>: Indicates if the solution of the previous condition should be used to initiate the next (only works when <code>condition=all</code>)
- <code>-v</code>, <code>--verbose</code>: Indicates if the program should run in verbose mode (can conflict with the option <code>-print</code>).
- <code>-vv</code>, <code>--extra-verbose</code>: Indicates if the program should run in extra-verbose mode (can conflict with the option <code>-print</code>).

# 5) SRC model format tutorial <a name="src_model_format_tutorial"></a>

A tutorial is available to better understand the content of a self-replicating cell model:

<p align="center">
<a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/src_model_format_tutorial.md" target="_blank">:link: SRC model format tutorial</a>
</p>

# 6) Units conversion tutorial <a name="units_conversion_tutorial"></a>

A tutorial is available for users starting from standard stoichiometric coefficients and kinetic parameters, and wanting to convert them into GBA formalism:

<p align="center">
<a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/units_conversion_tutorial.ipynb" target="_blank">:link: Units conversion tutorial</a>
</p>

# 7) Copyright <a name="copyright"></a>

Copyright © 2024-2025 Charles Rocabert.

# 8) License <a name="license"></a>

<p align="justify">
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
</p>

<p align="justify">
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
</p>

<p align="justify">
You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
</p>
