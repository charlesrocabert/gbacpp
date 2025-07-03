/**
 * \file      Macros.hpp
 * \author    Charles Rocabert
 * \date      12-08-2024
 * \copyright GBAcpp. Copyright © 2024-2025 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Definition of macro-variables
 */

/************************************************************************
 * GBAcpp (Growth Balance Analysis for C++)
 * Copyright © 2024-2025 Charles Rocabert
 * Web: https://github.com/charlesrocabert/gbacpp
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
 ************************************************************************/

#ifndef __GBAcpp__Macros__
#define __GBAcpp__Macros__

#define DECREASING_DT_FACTOR 5.0    /*!< Factor dividing dt                                      */
#define INCREASING_DT_FACTOR 2.0    /*!< Factor multiplying dt                                   */
#define INCREASING_DT_COUNT  100    /*!< Number of constant dt iterations to increase it         */
#define MIN_DT               1e-100 /*!< Minimal dt value                                        */
#define EXPORT_DATA_COUNT    100    /*!< Data is exported at this period in number of iterations */
#define REGULATION_SIGMA     10.0   /*!< Metabolite Gaussian kernel function width               */


#endif /* defined(__GBAcpp__Macros__) */

