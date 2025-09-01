/**
 * \file      Macros.hpp
 * \author    Charles Rocabert
 * \date      12-08-2024
 * \copyright gbacpp. Copyright © 2024-2025 Charles Rocabert. All rights reserved
 * \license   This work is licensed under the terms of the MIT license
 * \brief     Definition of macro-variables
 */

/********************************************************************************
 * gbacpp (growth balance analysis for C++)
 * Web: https://github.com/charlesrocabert/gbacpp
 *
 * MIT License
 *
 * Copyright © 2024-2025 Charles Rocabert. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ********************************************************************************/

#ifndef __gbacpp__Macros__
#define __gbacpp__Macros__


#define DECREASING_DT_FACTOR 5.0    /*!< Timestep dividing factor                                  */
#define INCREASING_DT_FACTOR 2.0    /*!< Timestep multiplying factor                               */
#define INCREASING_DT_COUNT  100    /*!< Iterations with constant timestep required to increase it */
#define MIN_DT               1e-100 /*!< Minimal timestep value                                    */
#define EXPORT_DATA_COUNT    100    /*!< Timestep window for data export                           */
#define REGULATION_SIGMA     10.0   /*!< Metabolite Gaussian kernel function width                 */


#endif /* defined(__gbacpp__Macros__) */

