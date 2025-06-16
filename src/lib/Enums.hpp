/**
 * \file      Enums.hpp
 * \author    Charles Rocabert
 * \date      21-08-2024
 * \copyright GBAcpp. Copyright © 2024-2025 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Definition of enumerations
 */

/****************************************************************************
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
 ****************************************************************************/

#ifndef __GBAcpp__Enums__
#define __GBAcpp__Enums__


/**
 * \brief   Reaction type
 * \details --
 */
enum rtype
{
  IMM   = 0, /*!< Irreversible Michaelis-Menten                                */
  IMMI  = 1, /*!< Irreversible Michaelis-Menten with inhibition                */
  IMMA  = 2, /*!< Irreversible Michaelis-Menten with activation                */
  IMMIA = 3, /*!< Irreversible Michaelis-Menten with inhibition and activation */
  IMMR  = 4, /*!< Irreversible Michaelis-Menten with regulation                */
  RMM   = 5  /*!< Reversible Michaelis-Menten                                  */
};

/**
 * \brief   Reaction direction
 * \details --
 */
enum rdir
{
  FORWARD    = 0, /*!< Forward irreversible reaction  */
  BACKWARD   = 1, /*!< Backward irreversible reaction */
  REVERSIBLE = 2  /*!< Reversible reaction            */
};


#endif /* defined(__GBAcpp__Enums__) */

