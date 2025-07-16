#!/bin/bash

#*******************************************************************************
# gbacpp (growth balance analysis for C++)
# Web: https://github.com/charlesrocabert/gbacpp
#
# MIT License
#
# Copyright © 2024-2025 Charles Rocabert
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#*******************************************************************************

echo "\n********************************************************************************"
echo "* gbacpp (growth balance analysis for C++)"
echo "* Web: https://github.com/charlesrocabert/gbacpp"
echo "*"
echo "* MIT License"
echo "*"
echo "* Copyright © 2024-2025 Charles Rocabert"
echo "*"
echo "* Permission is hereby granted, free of charge, to any person obtaining a copy"
echo "* of this software and associated documentation files (the \"Software\"), to deal"
echo "* in the Software without restriction, including without limitation the rights"
echo "* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell"
echo "* copies of the Software, and to permit persons to whom the Software is"
echo "* furnished to do so, subject to the following conditions:"
echo "*"
echo "* The above copyright notice and this permission notice shall be included in all"
echo "* copies or substantial portions of the Software."
echo "*"
echo "* THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR"
echo "* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,"
echo "* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE"
echo "* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER"
echo "* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,"
echo "* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE"
echo "* SOFTWARE."
echo "********************************************************************************\n"

cd cmake
bash make_clean.sh
cmake -DCMAKE_BUILD_TYPE=Release ..
make install

