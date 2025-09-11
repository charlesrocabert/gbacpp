#!/usr/bin/env bash

#***********************************************************************
# gbacpp (growth balance analysis for C++)
# Web: https://github.com/charlesrocabert/gbacpp
# Copyright © 2024-2025 Charles Rocabert.
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#***********************************************************************

find_model_optimum -path ../../gbamodels/CSV -name A -condition all -output ./output -verbose -optimum
find_model_optimum -path ../../gbamodels/CSV -name B -condition all -output ./output -verbose -optimum
find_model_optimum -path ../../gbamodels/CSV -name C -condition all -output ./output -verbose -optimum
find_model_optimum -path ../../gbamodels/CSV -name D -condition all -output ./output -verbose -optimum
find_model_optimum -path ../../gbamodels/CSV -name EC12b -condition all -output ./output -verbose -optimum
find_model_optimum -path ../../gbamodels/CSV -name BaseModel -condition all -output ./output -verbose -optimum
find_model_optimum -path ../../gbamodels/CSV -name ExtendedModel -condition all -output ./output -verbose -optimum
find_model_optimum -path ../../gbamodels/CSV -name EC3rev -condition all -output ./output -verbose -optimum 

