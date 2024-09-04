/**
 * \file      MMSYN_tests.cpp
 * \author    Charles Rocabert
 * \date      29-08-2024
 * \copyright GBA_Evolution. Copyright © 2024 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     MMSYN_tests executable
 */

/****************************************************************************
 * GBA_Evolution (Evolutionary Algorithms for Growth Balance Analysis)
 * Copyright © 2024 Charles Rocabert
 * Web: https://github.com/charlesrocabert/GBA_Evolution_2
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "../cmake/Config.h"

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <sys/stat.h>
#include <assert.h>

#include "./lib/Macros.hpp"
#include "./lib/Enums.hpp"
#include "./lib/Structs.hpp"
#include "./lib/Model.hpp"


/**
 * \brief    main function
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \return   \e int
 */
int main(int argc, char const** argv)
{
  Model* model = new Model("./csv_models", "MMSYN_FCR_noRNA_noDNA_nodUTPase", GENOME_SCALE);
  model->load_model();
  model->initialize_variables();
  model->compute_gradient_ascent_trajectory("1", 0.1, 200, true);
  
  delete model;
  model = NULL;
  return EXIT_SUCCESS;
}

