/**
 * \file      Macros.hpp
 * \author    Charles Rocabert
 * \date      12-08-2024
 * \copyright gbacpp. Copyright © 2024-2025 Charles Rocabert.
 * \license   GNU General Public License v3 (GPLv3)
 * \brief     Definition of macro-variables
 */

/************************************************************************
 * gbacpp (growth balance analysis for C++)
 * Web: https://github.com/charlesrocabert/gbacpp
 * Copyright © 2024-2025 Charles Rocabert.
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ************************************************************************/

#ifndef __gbacpp__Macros__
#define __gbacpp__Macros__


#define DECREASING_DT_FACTOR 5.0    /*!< Timestep dividing factor                                  */
#define INCREASING_DT_FACTOR 2.0    /*!< Timestep multiplying factor                               */
#define INCREASING_DT_COUNT  100    /*!< Iterations with constant timestep required to increase it */
#define MIN_DT               1e-100 /*!< Minimal timestep value                                    */
#define H_MIN                1e-8   /*!< Minimal absolute Hessian value                            */
#define KAPPA_MAX            1e+4   /*!< Maximal absolute inverse Hessian value                    */
#define EXPORT_DATA_COUNT    100    /*!< Timestep window for data export                           */


#endif /* defined(__gbacpp__Macros__) */

