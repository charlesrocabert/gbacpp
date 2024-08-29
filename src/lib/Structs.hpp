/**
 * \file      Structs.hpp
 * \author    Charles Rocabert
 * \date      28-08-2024
 * \copyright GBA_Evolution. Copyright © 2024 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Definition of structures
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

#ifndef __GBA_Evolution__Structs__
#define __GBA_Evolution__Structs__

#include <iostream>
#include <vector>
#include <cstring>

#include "Macros.hpp"
#include "Enums.hpp"


/**
 * \brief   Genetic unit struct
 * \details Defines the structure of a genetic unit
 */
typedef struct
{
  double LB; /* Lower boundary */
  double UB; /* Upper boundary */
} boundaries;


#endif /* defined(__GBA_Evolution__Structs__) */
