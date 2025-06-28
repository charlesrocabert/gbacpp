#!/usr/bin/env bash

#***********************************************************************
# GBAcpp (Growth Balance Analysis for C++)
# Copyright © 2024-2025 Charles Rocabert
# Web: https://github.com/charlesrocabert/GBAcpp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

mkdir GBAcpp
mkdir ./GBAcpp/output
cp -r ../build ./GBAcpp/.
cp -r ../cmake ./GBAcpp/.
cp -r ../src ./GBAcpp/.
cp -r ../csv_models ./GBAcpp/.
cp -r ../CMakeLists.txt ./GBAcpp/.
cp -r run_GBAcpp_jobs.py ./GBAcpp/.
zip -r GBAcpp.zip GBAcpp
scp GBAcpp.zip dam82xot@storage.hpc.rz.uni-duesseldorf.de:/gpfs/project/dam82xot/.
scp run_GBAcpp_jobs.py dam82xot@storage.hpc.rz.uni-duesseldorf.de:./
rm GBAcpp.zip
rm -rf GBAcpp
