/**
 * \file      Macros.hpp
 * \author    Charles Rocabert
 * \date      12-08-2024
 * \copyright GBA_Evolution. Copyright © 2024 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Definition of macro-variables
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

#ifndef __GBA_Evolution__Macros__
#define __GBA_Evolution__Macros__


#define FLUX_BOUNDARY              10.0  /*!< Flux absolute boundary value                                         */
#define MIN_PARAMETER              1e-10 /*!< Tolerance to test equality of parameters with zero                   */
#define MIN_CONCENTRATION          1e-10 /*!< Minimum concentration value                                          */
#define MIN_FLUX_FRACTION          1e-10 /*!< Minimum flux fraction                                                */
#define DENSITY_CONSTRAINT_TOL     1e-10  /*!< Density tolerance threshold (|1-rho| < tol)                          */
#define NEGATIVE_C_TOL             1e-10 /*!< Negative C tolerance threshold (C > -tol)                            */
#define NEGATIVE_P_TOL             1e-10 /*!< Negative P tolerance threshold (P > -tol)                            */
#define DECREASING_DT_FACTOR       5.0   /*!< Factor dividing dt                                                   */
#define INCREASING_DT_FACTOR       2.0   /*!< Factor multiplying dt                                                */
#define INCREASING_DT_COUNT        100   /*!< Number of constant dt iterations to increase it                       */
#define TRAJECTORY_STABLE_MU_COUNT 2000  /*!< Number of stable mu values required to consider a trajectory stable  */
#define TRAJECTORY_CONVERGENCE_TOL 1e-12  /*!< Analytical trajectory convergence tolerance                          */
#define EXPORT_DATA_COUNT          500  /*!< Data is exported at this period in number of iterations              */
#define MCMC_CONVERGENCE_TOL       1e-5  /*!< MCMC trajectory convergence tolerance                                */
#define POPLEVEL_CONVERGENCE_TOL   1e-5  /*!< Population-level trajectory convergence tolerance                    */
#define EFM_TOL                    1e-5  /*!< Tolerance threshold below which EFM values are considered to be zero */


#endif /* defined(__GBA_Evolution__Macros__) */
