<p align="center">
  <img src="https://github.com/user-attachments/assets/b60abce0-b7a2-4790-a3dc-7d059c8f37da" width=150 />
</p>
<h3 align="center">Growth Balance Analysis for C++</h3>

<p align="center">
<br />
<a href="https://github.com/charlesrocabert/gbacpp/releases/latest"><img src="https://img.shields.io/github/release/charlesrocabert/gbacpp/all.svg" /></a>&nbsp;
<a href="https://github.com/charlesrocabert/gbacpp/actions/workflows/cmake-multi-platform.yml"><img src="https://github.com/charlesrocabert/gbacpp/actions/workflows/cmake-multi-platform.yml/badge.svg" /></a>&nbsp;
<a href="https://github.com/charlesrocabert/gbacpp/LICENSE.html"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" /></a>
</p>

-----------------

# Roadmap

- [x] Gradient ascent (best for full column-rank CGMs<sup>*</sup>),
- [ ] MCMC algorithm,
- [ ] Forward-in-time population level simulations,
- [ ] Lineage tracking. 

<sup>*</sup>Cell growth model
<sup>*</sup>Markov chain Monte Carlo

# Installation instructions <a name="installation_instructions"></a>
Download the <a href="https://github.com/charlesrocabert/gbacpp/releases/latest">latest release</a> of <strong>gbacpp</strong>, and save it to a directory of your choice. Open a terminal and use the <code>cd</code> command to navigate to this directory. Then follow the steps below to compile and build the executables.

## Supported platforms <a name="supported_platforms"></a>
<strong>gbacpp</strong> software has been successfully tested on Unix/Linux and macOS platforms.

## Dependencies <a name="dependencies"></a>
* A C++ compiler (GCC, LLVM, ...),
* CMake (command line version),
* GSL for C/C++,
* CBLAS for C/C++.

## Software compilation <a name="software_compilation"></a>

### User mode
To compile <strong>gbacpp</strong>, run the following instructions on the command line:

    cd cmake/

and

    bash make.sh

### Debug mode
To compile the software in DEBUG mode, use <code>make_debug.sh</code> script instead of <code>make.sh</code>:

    bash make_debug.sh

This mode should only be used for test or development phases.

### Executable files emplacement
Binary executable files are in <code>build/bin</code> folder.

### Delete compiled files
To delete compiled files and binary executables, run:

    bash make_clean.sh
    
# First usage <a name="first_usage"></a>
Once <strong>gbacpp</strong> has been compiled, follow the next steps for a first usage of the software.

TO DO.

## Copyright <a name="copyright"></a>
Copyright © 2024-2025 Charles Rocabert. All rights reserved.

## License <a name="license"></a>
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
