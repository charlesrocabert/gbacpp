/**
 * \file      Model.hpp
 * \author    Charles Rocabert
 * \date      22-07-2024
 * \copyright gbacpp. Copyright © 2024-2025 Charles Rocabert. All rights reserved
 * \license   This work is licensed under the terms of the MIT license
 * \brief     Model class declaration
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


#ifndef __gbacpp__Model__
#define __gbacpp__Model__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>

#include "Macros.hpp"
#include "Enums.hpp"


class Model
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Model( void ) = delete;
  Model( std::string model_path, std::string model_name );
  Model( const Model& model ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Model( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  inline double get_mu( void );
  inline double get_doubling_time( void );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  Model& operator=(const Model&) = delete;
  
  inline void set_tol( double tol );
  inline void set_condition( std::string condition );
  inline void initialize_f( void );
  inline void calculate_f_from_f_trunc( void );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  
  void read_from_csv( void );
  void read_random_solutions( void );
  
  void compute_optimum( std::string condition, bool print_optimum, bool write_optimum, bool write_trajectory, std::string output_path, int stable_count, int max_iter, bool verbose );
  void compute_optimum_by_condition( bool print_optimum, bool write_optimum, bool write_trajectory, std::string output_path, int stable_count, int max_iter, bool verbose );
  //void compute_optimum_by_random_solution( std::string condition, bool print_optimum, bool write_optimum, bool write_trajectory, std::string output_path, int stable_count, int max_iter, bool verbose );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  
  bool is_path_exist( std::string path );
  bool is_file_exist( std::string filename );
  
  bool compute_gradient_ascent( std::string condition, bool write_trajectory, std::string output_path, int stable_count, int max_iter, bool verbose );
  
  void open_trajectory_output_files( std::string output_path, std::string condition );
  void write_trajectory_output_files( std::string condition, int iter, double t, double dt );
  void close_trajectory_ouput_files( void );
  
  void open_optimum_output_files( std::string output_path, std::string condition );
  void write_optimum_output_files( std::string condition, bool converged, double runtime );
  void close_optimum_ouput_files( void );
  
  void print_to_standard_ouput( std::string condition, bool converged, double runtime );
  
  void load_metabolite_identifiers( void );
  void load_reaction_identifiers( void );
  void load_vector_sizes( void );
  void load_M( void );
  void load_K( void );
  void load_KI( void );
  void load_KA( void );
  void load_KR( void );
  void load_kcat( void );
  void load_conditions( void );
  void load_constant_reactions( void );
  void load_f0( void );
  void reload_f0( void );
  
  void initialize_variables( void );
  void initialize_static_variables( void );
  void initialize_dynamic_variables( void );
  
  double gaussian_term( double x, double mu );
  double gaussian_kernel( double x, double mu );
  void calculate( void );
  void compute_c( void );
  void compute_xc( void );
  void iMM( int j );
  void iMMi( int j );
  void iMMa( int j );
  void iMMia( int j );
  void iMMr( int j );
  void rMM( int j );
  void compute_tau( int j );
  void diMM( int j );
  void diMMi( int j );
  void diMMa( int j );
  void diMMia( int j );
  void drMM( int j );
  void diMMr( int j );
  void compute_dtau( int j );
  void compute_mu( void );
  void compute_v( void );
  void compute_p( void );
  void compute_b( void );
  void compute_density( void );
  void compute_dmu_f( void );
  void compute_GCC_f( void );
  void calculate_first_order_terms( void );
  void calculate_second_order_terms( void );
  void check_model_consistency( void );
  void block_reactions( void );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*----------------------------------------------- Model path and name */
  
  std::string _model_path; /*!< Model path */
  std::string _model_name; /*!< Model name */
  
  /*----------------------------------------------- Tolerance value */
  
  double _tol; /*!< Tolerance value */
  
  /*----------------------------------------------- Identifier lists */
  
  std::vector<std::string>             _metabolite_ids;     /*!< List of all metabolite ids      */
  std::vector<std::string>             _x_ids;              /*!< List of external metabolite ids */
  std::vector<std::string>             _c_ids;              /*!< List of internal metabolite ids */
  std::vector<std::string>             _reaction_ids;       /*!< List of reaction ids            */
  std::vector<std::string>             _condition_ids;      /*!< List of condition ids           */
  std::vector<std::string>             _condition_params;   /*!< List of condition parameter ids */
  std::unordered_map<std::string, int> _metabolite_indices; /*!< Map of metabolite indices       */
  std::unordered_map<std::string, int> _reaction_indices;   /*!< Map of reaction indices         */
  std::unordered_map<std::string, int> _condition_indices;  /*!< Map of condition indices        */
  
  /*----------------------------------------------- Model structure */
  
  gsl_matrix*                             _Mx;                 /*!< Total mass fraction matrix              */
  gsl_matrix*                             _M;                  /*!< Internal mass fraction matrix           */
  gsl_matrix*                             _K;                  /*!< Complete K matrix                       */
  gsl_matrix*                             _KM_f;               /*!< Forward KM matrix                       */
  gsl_matrix*                             _KM_b;               /*!< Backward KM matrix                      */
  gsl_matrix*                             _KI;                 /*!< KI matrix                               */
  gsl_matrix*                             _KA;                 /*!< KA matrix                               */
  gsl_matrix*                             _KR;                 /*!< KR matrix                               */
  gsl_vector*                             _kcat_f;             /*!< Forward kcat vector                     */
  gsl_vector*                             _kcat_b;             /*!< Backward kcat vector                    */
  rtype*                                  _type;               /*!< Reaction type                           */
  gsl_matrix*                             _conditions;         /*!< List of conditions                      */
  std::unordered_map<std::string, double> _constant_reactions; /*!< Constant reactions with constant values */

  /*----------------------------------------------- Vector lengths */
  
  int _nx; /*!< Number of external metabolites */
  int _nc; /*!< Number of internal metabolites */
  int _ni; /*!< Total number of metabolites    */
  int _nj; /*!< Number of reactions            */
  
  /*----------------------------------------------- Other model variables */
  
  gsl_vector* _sM; /*!< Columns sum of M                   */
  int         _r;  /*!< Ribosome reaction index            */
  int         _a;  /*!< Total proteins concentration index */
  
  /*----------------------------------------------- GBA external condition variables */
  
  std::string _condition; /*!< Current environmental condition    */
  double      _rho;       /*!< Current total density              */
  gsl_vector* _x;         /*!< External metabolite concentrations */
  
  /*----------------------------------------------- GBA first order variables */
  
  gsl_vector* _f0;                    /*!< Initial state                      */
  gsl_vector* _f;                     /*!< Flux fractions vector              */
  gsl_vector* _f_trunc;               /*!< Truncated flux fractions vector    */
  gsl_vector* _c;                     /*!< Internal metabolite concentrations */
  gsl_vector* _xc;                    /*!< Metabolite concentrations          */
  gsl_vector* _tau_j;                 /*!< Tau values (turnover times)        */
  gsl_vector* _v;                     /*!< Fluxes vector                      */
  gsl_vector* _p;                     /*!< Protein concentrations vector      */
  gsl_vector* _b;                     /*!< Biomass fractions vector           */
  double      _density;               /*!< Cell's relative density            */
  double      _mu;                    /*!< Growth rate                        */
  double      _doubling_time;         /*!< Doubling time                      */
  bool        _consistent;            /*!< Is the model consistent?           */
  bool        _adjust_concentrations; /*!< Adjust concentration vector c      */
  
  /*----------------------------------------------- GBA second order variables */
  
  gsl_matrix* _ditau_j; /*!< Tau derivative values                               */
  gsl_vector* _dmu_f;   /*!< Local mu derivatives with respect to f              */
  gsl_vector* _GCC_f;   /*!< Local growth control coefficients with respect to f */
  
  /*----------------------------------------------- Variables for calculation optimization */
  
  gsl_vector_view _x_view;      /*!< x segment view of vector xc           */
  gsl_vector_view _c_view;      /*!< c segment view of vector xc           */
  double          _dmu_f_term1; /*!< Variable for the calculation of dmu_f */
  gsl_vector*     _dmu_f_term2; /*!< Variable for the calculation of dmu_f */
  gsl_matrix*     _dmu_f_term3; /*!< Variable for the calculation of dmu_f */
  gsl_vector*     _dmu_f_term4; /*!< Variable for the calculation of dmu_f */
  gsl_vector*     _dmu_f_term5; /*!< Variable for the calculation of dmu_f */
  double          _mu_diff;     /*!< Next mu to current mu differential    */
  
  /*----------------------------------------------- Solutions */
  
  int                                  _nb_random_solutions; /*!< Number of random solutions                    */
  std::unordered_map<int, gsl_vector*> _random_solutions;    /*!< List of random f vectors                      */
  
  /*----------------------------------------------- Output files */
  
  std::ofstream _state_trajectory_file; /*!< Model trajectory output file    */
  std::ofstream _f_trajectory_file;     /*!< f vector trajectory output file */
  std::ofstream _c_trajectory_file;     /*!< c vector trajectory output file */
  std::ofstream _v_trajectory_file;     /*!< v vector trajectory output file */
  std::ofstream _p_trajectory_file;     /*!< p vector trajectory output file */
  std::ofstream _b_trajectory_file;     /*!< b vector trajectory output file */
  
  std::ofstream _state_optimum_file; /*!< Model optimum output file    */
  std::ofstream _f_optimum_file;     /*!< f vector optimum output file */
  std::ofstream _c_optimum_file;     /*!< c vector optimum output file */
  std::ofstream _v_optimum_file;     /*!< v vector optimum output file */
  std::ofstream _p_optimum_file;     /*!< p vector optimum output file */
  std::ofstream _b_optimum_file;     /*!< b vector optimum output file */
};

/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get mu
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Model::get_mu( void )
{
  return _mu;
}

/**
 * \brief    Get doubling time
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Model::get_doubling_time( void )
{
  return _doubling_time;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/**
 * \brief    Set the tolerance value
 * \details  --
 * \param    double tol
 * \return   \e void
 */
inline void Model::set_tol( double tol )
{
  if (tol <= 0.0)
  {
    throw std::invalid_argument("> Error: tolerance value is too low ("+std::to_string(tol)+")");
  }
  _tol = tol;
}

/**
 * \brief    Set the external condition
 * \details  --
 * \param    std::string condition
 * \return   \e void
 */
inline void Model::set_condition( std::string condition )
{
  if (_condition_indices.find(condition) == _condition_indices.end())
  {
    throw std::invalid_argument("> Error: Unknown condition "+condition);
  }
  _condition = condition;
  int cond_index     = _condition_indices[condition];
  for(int i = 0; i < (int)_condition_params.size(); i++)
  {
    std::string param = _condition_params[i];
    if (param == "rho")
    {
      _rho = gsl_matrix_get(_conditions, i, cond_index);
    }
    else
    {
      assert(_metabolite_indices.find(param) != _metabolite_indices.end());
      int met_index = _metabolite_indices[param];
      gsl_vector_set(_x, met_index, gsl_matrix_get(_conditions, i, cond_index));
    }
  }
}

/**
 * \brief    Initialize the vectors f and f_trunc from f0
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Model::initialize_f( void )
{
  gsl_vector_memcpy(_f, _f0);
  gsl_vector_view f_view = gsl_vector_subvector(_f, 1, _nj-1);
  gsl_vector_memcpy(_f_trunc, &f_view.vector);
}

/**
 * \brief    Calculate f from f_trunc
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Model::calculate_f_from_f_trunc( void )
{
  gsl_vector_view sM_view = gsl_vector_subvector(_sM, 1, _nj-1);
  double term = 0.0;
  gsl_blas_ddot(&sM_view.vector, _f_trunc, &term);
  term = (1.0-term)/gsl_vector_get(_sM, 0);
  gsl_vector_view f_view = gsl_vector_subvector(_f, 1, _nj-1);
  gsl_vector_memcpy(&f_view.vector, _f_trunc);
  gsl_vector_set(_f, 0, term);
}


#endif /* defined(__gbacpp__Model__) */
