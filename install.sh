#!/bin/bash

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

echo "\n************************************************************************"
echo "* gbacpp (growth balance analysis for C++)"
echo "* Web: https://github.com/charlesrocabert/gbacpp"
echo "* Copyright © 2024-2025 Charles Rocabert."
echo "*"
echo "* This program is free software: you can redistribute it and/or modify"
echo "* it under the terms of the GNU General Public License as published by"
echo "* the Free Software Foundation, either version 3 of the License, or"
echo "* (at your option) any later version."
echo "*"
echo "* This program is distributed in the hope that it will be useful,"
echo "* but WITHOUT ANY WARRANTY; without even the implied warranty of"
echo "* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
echo "* GNU General Public License for more details."
echo "*"
echo "* You should have received a copy of the GNU General Public License"
echo "* along with this program.  If not, see <https://www.gnu.org/licenses/>."
echo "************************************************************************\n"

cd cmake
bash make_clean.sh
cmake -DCMAKE_BUILD_TYPE=Release ..
make install

