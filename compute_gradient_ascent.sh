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

./build/bin/compute_gradient_ascent -path /Users/charlesrocabert/git/charlesrocabert/gbapy/tutorials/MMSYN_tutorial/models \
-name mmsyn_fcr_v2 -condition 67 -dt 0.01 -maxt 100000 -max-mu-count 1000 -save -output /Users/charlesrocabert/git/charlesrocabert/gbapy/tutorials/MMSYN_tutorial/output/optimization

