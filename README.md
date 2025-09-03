<p align="center">
  <img src="https://github.com/user-attachments/assets/3ba498cc-a202-49b1-ade1-ac0488dd30b9" width=200 />

</p>
<h3 align="center">Growth Balance Analysis for C++</h3>

<p align="center">
<br />
<a href="https://github.com/charlesrocabert/gbacpp/releases/latest"><img src="https://img.shields.io/github/release/charlesrocabert/gbacpp/all.svg" /></a>&nbsp;
<a href="https://github.com/charlesrocabert/gbacpp/actions/workflows/cmake-multi-platform.yml"><img src="https://github.com/charlesrocabert/gbacpp/actions/workflows/cmake-multi-platform.yml/badge.svg" /></a>&nbsp;
<a href="https://github.com/charlesrocabert/gbacpp/LICENSE.html"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" /></a>
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
The optimization process relies on a gradient ascent approach, and is preferred for models offering a convex solution space (typically, when all reactions in the metabolic network are linearly independent).
</p>

<p align="justify">
:ballot_box_with_check: Note that SRC models must comply to a standardized format. Guidelines are available in the <a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/cgm_format_tutorial.md" target="_blank">SRC model format tutorial</a>.

:ballot_box_with_check: When building a SRC model, stoichiometric coefficients, and kinetic parameters must be converted following GBA formalism. See the <a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/units_conversion_tutorial.ipynb" target="_blank">units conversion tutorial</a>.

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
  - [4.1) Find an optimum](#find_optimum)
  - [4.2) Optimization parameters](#optimization_parameters)
  - [4.3) Usage example](#usage_example)
  - [4.4) Ready-to-use examples](#examples)
- [5) SRC model format tutorial](#src_model_format_tutorial)
- [6) Units conversion tutorial](#units_conversion_tutorial)
- [7) Copyright](#copyright)
- [8) License](#license)

# 1) Roadmap <a name="roadmap"></a>

| Task | Status |
|---|---|
| Gradient ascent (best for full column-rank SRC models) | <a href="https://postgresql.org"><img src="https://img.shields.io/badge/Status-Done-green.svg"/></a> |
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

### 4.1) Find an optimum <a name="find_optimum"></a>
To run a gradient ascent optimization on a SRC model, execute the following command line:

    find_model_optimum <parameters>

The command line parameters are described below. The description is also available by executing the following command line in a terminal:

    find_model_optimum -h

### 4.2) Optimization parameters <a name="optimization_parameters"></a>

- <code>-h</code>, <code>--help</code>: Print the help, then exit,
- <code>-v</code>, <code>--version</code>: Print the current version, then exit,

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
- <code>-stable</code>, <code>--stable-count</code>: Specify the maximal number of iterations with unchanged growth rate as a stop criterium ($10,000$ by default),
- <code>-max</code>, <code>--max-iter</code>: Specify the maximal number of iterations as a stop criterium ($100,000,000$ by default),
- <code>-verbose</code>, <code>--verbose</code>: Indicates if the program should run in verbose mode (can conflict with the option <code>-print</code>).
        
### 4.3) Usage example <a name="usage_example"></a>

In this example, growth rate optimums are calculated for all the external conditions of the <em>Escherichia coli</em> toy model EC12b (see folder <code>./examples/toy_models</code>).

First, navigate to the folder <code>./examples</code> using the <code>cd</code> command:

    cd ./examples

Then, call the optimization algorithm:

    find_model_optimum -path ./toy_models -name EC12b -condition all -output ./output -verbose

Here, optimums are calculated for all conditions, and saved in the folder <code>./examples/output</code>. Verbose mode is activated to get insights in the optimization process.

### 4.4) Ready-to-use examples <a name="examples"></a>
Ready-to-use examples are available in the folder <code>./examples</code> (navigate to the folder <code>examples</code> using the <code>cd</code> command):

• <code>model_A_condition_1.sh</code>: This script will run a single gradient ascent on model A in external condition 1 (2 reactions, 2 metabolites). You can execute it using the following command line:

    bash model_A_condition_1.sh

At the end of the optimization, CSV files are written in the folder <code>./examples/output</code>. You can edit the parameter values at will to test the behaviour of the gradient ascent. See below for a full description of the parameters.

• <code>model_B_all_conditions.sh</code>: This script will calculate the optimal growth rate on model B for all external conditions. You can execute it using the following command line:

    bash model_B_all_conditions.sh

All the optimums are written in the folder <code>./examples/output</code>.

## 5) SRC model format tutorial <a name="src_model_format_tutorial"></a>

A tutorial is available to better understand the content of cell growth model:

<p align="center">
<a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/src_model_format_tutorial.md" target="_blank">:link: SRC model format tutorial</a>
</p>

## 6) Units conversion tutorial <a name="units_conversion_tutorial"></a>

A tutorial is available for users starting from standard metabolic models stoichiometric coefficients and kinetic parameters, and wanting to convert them into GBA formalism:

<p align="center">
<a href="https://github.com/charlesrocabert/gbacpp/blob/main/tutorials/units_conversion_tutorial.ipynb" target="_blank">:link: Units conversion tutorial</a>
</p>

## 7) Copyright <a name="copyright"></a>
Copyright © 2024-2025 Charles Rocabert. All rights reserved.

## 8) License <a name="license"></a>

MIT License

Copyright (c) 2024-2025 Charles Rocabert

<p align="justify">
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
</p>

<p align="justify">
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
</p>

<p align="justify">
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
</p>
