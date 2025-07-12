#!/usr/bin/env bash

#***********************************************************************
# gbacpp (growth balance analysis for C++)
# Copyright © 2024-2025 Charles Rocabert
# Web: https://github.com/charlesrocabert/gbacpp
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

find_cgm_optimum -path ./models -name A -condition all -output ./output -verbose
find_cgm_optimum -path ./models -name B -condition all -output ./output -verbose
find_cgm_optimum -path ./models -name C -condition all -output ./output -verbose
find_cgm_optimum -path ./models -name D -condition all -output ./output -verbose
find_cgm_optimum -path ./models -name EC12b -condition all -output ./output -verbose
find_cgm_optimum -path ./models -name BaseModel -condition all -output ./output -verbose
find_cgm_optimum -path ./models -name ExtendedModel -condition all -output ./output -verbose
find_cgm_optimum -path ./models -name EC3rev -condition all -output ./output -verbose

