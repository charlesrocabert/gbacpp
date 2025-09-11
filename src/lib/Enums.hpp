/**
 * \file      Enums.hpp
 * \author    Charles Rocabert
 * \date      21-08-2024
 * \copyright gbacpp. Copyright © 2024-2025 Charles Rocabert.
 * \license   GNU General Public License v3 (GPLv3)
 * \brief     Definition of enumerations
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

#ifndef __gbacpp__Enums__
#define __gbacpp__Enums__


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
  RMM   = 4  /*!< Reversible Michaelis-Menten                                  */
};


#endif /* defined(__gbacpp__Enums__) */

