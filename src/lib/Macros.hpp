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


#define MIN_PARAMETER                 1e-10  /*!< Tolerance to test equality of parameters with zero                  */
#define MIN_CONCENTRATION             1e-10  /*!< Minimum concentration value                                         */
#define MIN_FLUX_FRACTION             1e-10  /*!< Minimum flux fraction                                               */
#define DENSITY_CONSTRAINT_TOL        1e-10  /*!< Density tolerance threshold (|1-rho| < tol)                         */
#define NEGATIVE_C_TOL                1e-10  /*!< Negative C tolerance threshold (C > -tol)                           */
#define NEGATIVE_P_TOL                1e-10  /*!< Negative P tolerance threshold (P > -tol)                           */
#define DECREASING_DT_FACTOR          5.0    /*!< Factor dividing dt                                                  */
#define INCREASING_DT_FACTOR          2.0    /*!< Factor multiplying dt                                               */
#define INCREASING_DT_COUNT           100    /*!< Number of constant dt iterations to increase it                     */
#define MIN_DT                        1e-100 /*!< Minimal dt value                                                    */
#define TRAJECTORY_STABLE_MU_COUNT    10000  /*!< Number of stable mu values required to consider a trajectory stable */
#define TRAJECTORY_INCONSISTENT_COUNT 10000  /*!< Number of inconsistent trajectories required to stop the simulation */
#define TRAJECTORY_CONVERGENCE_TOL    1e-10  /*!< Analytical trajectory convergence tolerance                         */
#define EXPORT_DATA_COUNT             100    /*!< Data is exported at this period in number of iterations             */


#endif /* defined(__GBAcpp__Macros__) */

