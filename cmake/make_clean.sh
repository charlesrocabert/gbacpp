#!/bin/bash

#***************************************************************************
# GBA_Evolution (Evolutionary Algorithms for Growth Balance Analysis)
# Copyright © 2024 Charles Rocabert
# Web: https://github.com/charlesrocabert/GBA_Evolution_2
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
# ***************************************************************************

make clean
rm -rf ../build/bin/*
rm -rf CMakeFiles
rm cmake_install.cmake
rm CMakeCache.txt
rm Config.h
rm Makefile
