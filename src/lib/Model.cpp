/**
 * \file      Model.cpp
 * \author    Charles Rocabert
 * \date      22-07-2024
 * \copyright gbacpp. Copyright © 2024-2025 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Model class definition
 */

/****************************************************************************
 * gbacpp (growth balance analysis for C++)
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

#include "Model.hpp"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    std::string model_path
 * \param    std::string model_name
 * \return   \e void
 */
Model::Model( std::string model_path, std::string model_name )
{
  /*----------------------------------------------- Model path and name */
  
  assert(is_path_exist(model_path +"/"+model_name));
  _model_path = model_path;
  _model_name = model_name;
  
  /*----------------------------------------------- Tolerance value */
  
  _tol = 1e-10;
  
  /*----------------------------------------------- Identifier lists */
  
  _metabolite_ids.clear();
  _x_ids.clear();
  _c_ids.clear();
  _reaction_ids.clear();
  _condition_ids.clear();
  _condition_params.clear();
  _metabolite_indices.clear();
  _reaction_indices.clear();
  
  /*----------------------------------------------- Model structure */
  
  _Mx         = NULL;
  _M          = NULL;
  _KM         = NULL;
  _KM_f       = NULL;
  _KM_b       = NULL;
  _KI         = NULL;
  _KA         = NULL;
  _KR         = NULL;
  _kcat_f     = NULL;
  _kcat_b     = NULL;
  _type       = NULL;
  _conditions = NULL;
  _constant_reactions.clear();

  /*----------------------------------------------- Vector lengths */
  
  _nx = 0;
  _nc = 0;
  _ni = 0;
  _nj = 0;
  
  /*----------------------------------------------- Other model variables */
  
  _sM = NULL;
  _r  = 0;
  _a  = 0;
  
  /*----------------------------------------------- GBA external condition variables */
  
  _current_condition = "";
  _current_rho       = 0.0;
  _x                 = NULL;
  
  /*----------------------------------------------- GBA first order variables */
  
  _f0                    = NULL;
  _f                     = NULL;
  _f_trunc               = NULL;
  _c                     = NULL;
  _xc                    = NULL;
  _tau_j                 = NULL;
  _v                     = NULL;
  _p                     = NULL;
  _b                     = NULL;
  _density               = 0.0;
  _mu                    = 0.0;
  _doubling_time         = 0.0;
  _consistent            = false;
  _adjust_concentrations = false;
  
  /*----------------------------------------------- GBA second order variables */
  
  _ditau_j = NULL;
  _dmu_f   = NULL;
  _GCC_f   = NULL;
  
  /*----------------------------------------------- Variables for calculation optimization */
  
  _dmu_f_term1 = 0.0;
  _dmu_f_term2 = NULL;
  _dmu_f_term3 = NULL;
  _dmu_f_term4 = NULL;
  _dmu_f_term5 = NULL;
  
  /*----------------------------------------------- Solutions */
  
  _nb_random_solutions = 0;
  _random_solutions.clear();
  
  /*----------------------------------------------- Load the model */
  
  read_from_csv();
  initialize_variables();
}

/*----------------------------
 * DESTRUCTORS
 *----------------------------*/

/**
 * \brief    Destructor
 * \details  --
 * \param    void
 * \return   \e void
 */
Model::~Model( void )
{
  /*----------------------------------------------- Identifier lists */
  
  _metabolite_ids.clear();
  _x_ids.clear();
  _c_ids.clear();
  _reaction_ids.clear();
  _condition_ids.clear();
  _condition_params.clear();
  _metabolite_indices.clear();
  _reaction_indices.clear();
  
  /*----------------------------------------------- Model structure */
  
  gsl_matrix_free(_Mx);
  gsl_matrix_free(_M);
  gsl_matrix_free(_KM);
  gsl_matrix_free(_KM_f);
  gsl_matrix_free(_KM_b);
  gsl_matrix_free(_KI);
  gsl_matrix_free(_KA);
  gsl_matrix_free(_KR);
  gsl_vector_free(_kcat_f);
  gsl_vector_free(_kcat_b);
  delete[] _type;
  gsl_matrix_free(_conditions);
  _Mx         = NULL;
  _M          = NULL;
  _KM         = NULL;
  _KM_f       = NULL;
  _KM_b       = NULL;
  _KI         = NULL;
  _KA         = NULL;
  _KR         = NULL;
  _kcat_f     = NULL;
  _kcat_b     = NULL;
  _type       = NULL;
  _conditions = NULL;
  _constant_reactions.clear();
  
  /*----------------------------------------------- Other model variables */
  
  gsl_vector_free(_sM);
  _sM = NULL;
  
  /*----------------------------------------------- GBA external condition variables */
  
  gsl_vector_free(_x);
  _x = NULL;
  
  /*----------------------------------------------- GBA first order variables */
  
  gsl_vector_free(_f0);
  gsl_vector_free(_f);
  gsl_vector_free(_f_trunc);
  gsl_vector_free(_c);
  gsl_vector_free(_xc);
  gsl_vector_free(_tau_j);
  gsl_vector_free(_v);
  gsl_vector_free(_p);
  gsl_vector_free(_b);
  _f0         = NULL;
  _f          = NULL;
  _f_trunc    = NULL;
  _c          = NULL;
  _xc         = NULL;
  _tau_j      = NULL;
  _v          = NULL;
  _p          = NULL;
  _b          = NULL;
  
  /*----------------------------------------------- GBA second order variables */
  
  gsl_matrix_free(_ditau_j);
  gsl_vector_free(_dmu_f);
  gsl_vector_free(_GCC_f);
  _ditau_j = NULL;
  _dmu_f   = NULL;
  _GCC_f   = NULL;
  
  /*----------------------------------------------- Variables for calculation optimization */
  
  gsl_vector_free(_dmu_f_term2);
  gsl_matrix_free(_dmu_f_term3);
  gsl_vector_free(_dmu_f_term4);
  gsl_vector_free(_dmu_f_term5);
  _dmu_f_term2 = NULL;
  _dmu_f_term3 = NULL;
  _dmu_f_term4 = NULL;
  _dmu_f_term5 = NULL;
  
  /*----------------------------------------------- Solutions */
  
  for (int i = 0; i < _nb_random_solutions; i++)
  {
    gsl_vector_free(_random_solutions[i]);
    _random_solutions[i] = NULL;
  }
  _nb_random_solutions = 0;
  _random_solutions.clear();
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Read the model from CSV files
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::read_from_csv( void )
{
  load_metabolite_identifiers();
  load_reaction_identifiers();
  load_vector_sizes();
  load_M();
  load_KM();
  load_KI();
  load_KA();
  load_KR();
  load_kcat();
  load_conditions();
  load_constant_reactions();
  load_f0();
}

/**
 * \brief    Read pre-generated random solutions
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::read_random_solutions( void )
{
  assert(_random_solutions.size()==0);
  assert(is_path_exist(_model_path+"/"+_model_name));
  assert(is_file_exist(_model_path+"/"+_model_name+"/random_solutions.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/random_solutions.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string str_value;
  getline(file, line);
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Load the header and parse reaction indices */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  std::vector<std::string> header;
  std::stringstream flux(line.c_str());
  while(getline(flux, str_value, ';'))
  {
    header.push_back(str_value);
  }
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Parse each random solution                 */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  _nb_random_solutions = 0;
  while(getline(file, line))
  {
    /*** 2.1) Initialize the new solution vector ***/
    _random_solutions[_nb_random_solutions] = gsl_vector_alloc(_nj);
    gsl_vector_set_zero(_random_solutions[_nb_random_solutions]);
    /*** 2.2) Parse the line ***/
    std::stringstream flux(line.c_str());
    int pos = 0;
    while(getline(flux, str_value, ';'))
    {
      if (_reaction_indices.find(header[pos]) != _reaction_indices.end())
      {
        gsl_vector_set(_random_solutions[_nb_random_solutions], _reaction_indices[header[pos]], stod(str_value));
      }
      pos++;
    }
    /*** 2.3) Update the number of random solutions ***/
    _nb_random_solutions++;
  }
  file.close();
}

/**
 * \brief    Compute the optimum for one condition
 * \details  --
 * \param    std::string condition
 * \param    bool print_optimum
 * \param    bool write_trajectory
 * \param    std::string output_path
 * \param    int stable_count
 * \param    double max_t
 * \param    bool verbose
 * \return   \e bool
 */
void Model::compute_optimum( std::string condition, bool print_optimum, bool write_trajectory, std::string output_path, int stable_count, double max_t, bool verbose )
{
  std::clock_t begin = clock();
  bool converged     = compute_gradient_ascent(condition, write_trajectory, output_path, stable_count, max_t, verbose);
  std::clock_t end   = clock();
  double runtime     = double(end-begin)/CLOCKS_PER_SEC;
  open_optimum_output_files(output_path, condition);
  write_optimum_output_files(condition, converged, runtime);
  close_optimum_ouput_files();
  if (print_optimum)
  {
    print_to_standard_ouput(condition, converged, runtime);
  }
  if (verbose && converged)
  {
    std::cout << "> Condition " << condition << ": convergence reached (mu=" << _mu << ", runtime=" << runtime << ")" << std::endl;
  }
  else if (verbose && !converged)
  {
    std::cout << "> Condition " << condition << ": convergence not reached after T=" << max_t << " (mu=" << _mu << ", runtime=" << runtime << ")" << std::endl;
  }
}

/**
 * \brief    Compute the optimum for all the conditions
 * \details  --
 * \param    bool print_optimum
 * \param    bool write_trajectory
 * \param    std::string output_path
 * \param    int stable_count
 * \param    double max_t
 * \param    bool verbose
 * \return   \e bool
 */
void Model::compute_optimum_by_condition( bool print_optimum, bool write_trajectory, std::string output_path, int stable_count, double max_t, bool verbose )
{
  open_optimum_output_files(output_path, "all");
  for (int i = 0; i < (int)_condition_ids.size(); i++)
  {
    std::clock_t begin     = clock();
    std::string  condition = _condition_ids[i];
    bool         converged = compute_gradient_ascent(condition, write_trajectory, output_path, stable_count, max_t, verbose);
    std::clock_t end       = clock();
    double       runtime   = double(end-begin)/CLOCKS_PER_SEC;
    write_optimum_output_files(condition, converged, runtime);
    if (print_optimum)
    {
      print_to_standard_ouput(condition, converged, runtime);
    }
    if (verbose && converged)
    {
      std::cout << "> Condition " << condition << ": convergence reached (mu=" << _mu << ", runtime=" << runtime << ")" << std::endl;
    }
    else if (verbose && !converged)
    {
      std::cout << "> Condition " << condition << ": convergence not reached after T=" << max_t << " (mu=" << _mu << ", runtime=" << runtime << ")" << std::endl;
    }
  }
  close_optimum_ouput_files();
}

/**
 * \brief    Compute the optimum for all random solution and for a given condition
 * \details  --
 * \param    std::string condition
 * \param    bool print_optimum
 * \param    bool write_trajectory
 * \param    std::string output_path
 * \param    int stable_count
 * \param    double max_t
 * \param    bool verbose
 * \return   \e bool
 */
/*
void Model::compute_optimum_by_random_solution( std::string condition, bool print_optimum, bool write_trajectory, std::string output_path, int stable_count, double max_t, bool verbose )
{
  open_optimum_output_files(output_path, "random");
  for (int i = 0; i < _nb_random_solutions; i++)
  {
    
    std::clock_t begin = clock();
    gsl_vector_memcpy(_f0, _random_solutions[i]);
    bool         converged = compute_gradient_ascent(condition, write_trajectory, output_path, stable_count, max_t, verbose);
    std::clock_t end       = clock();
    double       runtime   = double(end-begin)/CLOCKS_PER_SEC;
    write_optimum_output_files(condition, converged, runtime);
    if (print_optimum)
    {
      print_to_standard_ouput(condition, converged, runtime);
    }
    if (verbose)
    {
      std::cout << "> Elapsed time for random solution " << i << ": " << runtime << " seconds" << std::endl;
    }
  }
  close_optimum_ouput_files();
}
*/

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Test path existence
 * \details  --
 * \param    std::string path
 * \return   \e bool
 */
bool Model::is_path_exist( std::string path )
{
  return std::filesystem::exists(path);
}

/**
 * \brief    Test file existence
 * \details  --
 * \param    std::string filename
 * \return   \e bool
 */
bool Model::is_file_exist( std::string filename )
{
  bool is_file_exist = false;
  std::ifstream infile(filename.c_str());
  is_file_exist = infile.good();
  infile.close();
  return is_file_exist;
}

/**
 * \brief    Compute the gradient ascent trajectory
 * \details  --
 * \param    std::string condition
 * \param    bool write_trajectory
 * \param    std::string output_path
 * \param    int stable_count
 * \param    double max_t
 * \param    bool verbose
 * \return   \e bool
 */
bool Model::compute_gradient_ascent( std::string condition, bool write_trajectory, std::string output_path, int stable_count, double max_t, bool verbose )
{
  auto it = std::find(_condition_ids.begin(), _condition_ids.end(), condition);
  if (it==_condition_ids.end())
  {
    throw std::invalid_argument("> Error: Unknown condition");
  }
  if (!is_path_exist(output_path))
  {
    throw std::invalid_argument("> Error: Path "+output_path+" does not exist");
  }
  if (stable_count < 0)
  {
    throw std::invalid_argument("> Error: The stable count parameter must be positive or null");
  }
  if (max_t <= 0.0)
  {
    throw std::invalid_argument("> Error: The maximum time must be positive");
  }
  if (write_trajectory)
  {
    open_trajectory_output_files(output_path, condition);
  }
  _adjust_concentrations = false;
  set_condition(condition);
  initialize_f();
  calculate();
  /*
  if (verbose)
  {
    std::cout << "> Initial growth rate = " << _mu << "\n";
  }
   */
  if (!_consistent)
  {
    throw std::runtime_error("> Error: The initial solution f0 is not consistent");
  }
  gsl_vector* previous_f_trunc = gsl_vector_alloc(_nj-1);
  gsl_vector* scaled_dmudt     = gsl_vector_alloc(_nj-1);
  gsl_vector_memcpy(previous_f_trunc, _f_trunc);
  double previous_mu         = 0.0;
  double t                   = 0.0;
  double dt                  = 0.01;
  int    dt_counter          = 0;
  int    constant_mu_counter = 0;
  int    nb_iterations       = 0;
  int    nb_successes        = 0;
  int    inconsistent_count  = 0;
  write_trajectory_output_files(condition, nb_iterations, t, dt);
  while (t < max_t)
  {
    nb_iterations++;
    if (constant_mu_counter >= stable_count)
    {
      break;
    }
    previous_mu = _mu;
    block_reactions();
    gsl_vector_view dmudt = gsl_vector_subvector(_GCC_f, 1, _nj-1);
    gsl_vector_memcpy(scaled_dmudt, &dmudt.vector);
    gsl_vector_scale(scaled_dmudt, dt);
    gsl_vector_add(_f_trunc, scaled_dmudt);
    calculate_f_from_f_trunc();
    calculate();
    if (_consistent && _mu >= previous_mu)
    {
      gsl_vector_memcpy(previous_f_trunc, _f_trunc);
      nb_successes++;
      inconsistent_count = 0;
      t                  = t+dt;
      dt_counter++;
      if (write_trajectory && nb_iterations%EXPORT_DATA_COUNT == 0)
      {
        if (verbose)
        {
          std::cout << " > " << nb_iterations << " iterations, " << nb_successes << " successes (mu=" << _mu << ", constant mu iters=" << constant_mu_counter << ", dt=" << dt << ")" << std::endl;
        }
        write_trajectory_output_files(condition, nb_iterations, t, dt);
      }
      if (fabs(_mu-previous_mu) < stable_count)
      {
        constant_mu_counter++;
      }
      else
      {
        constant_mu_counter--;
        if (constant_mu_counter < 0)
        {
          constant_mu_counter = 0;
        }
      }
      if (dt_counter == INCREASING_DT_COUNT)
      {
        dt         *= INCREASING_DT_FACTOR;
        dt_counter  = 0;
      }
    }
    else
    {
      gsl_vector_memcpy(_f_trunc, previous_f_trunc);
      calculate_f_from_f_trunc();
      calculate();
      assert(_consistent);
      inconsistent_count++;
      dt         /= DECREASING_DT_FACTOR;
      dt_counter  = 0;
    }
  }
  gsl_vector_free(previous_f_trunc);
  gsl_vector_free(scaled_dmudt);
  previous_f_trunc = NULL;
  scaled_dmudt     = NULL;
  if (write_trajectory)
  {
    write_trajectory_output_files(condition, nb_iterations, t, dt);
    close_trajectory_ouput_files();
  }
  if (constant_mu_counter >= stable_count)
  {
    return(true);
  }
  else
  {
    return(false);
  }
}

/**
 * \brief    Open trajectory output files
 * \details  Also writes headers
 * \param    std::string output_path
 * \param    std::string condition
 * \return   \e void
 */
void Model::open_trajectory_output_files( std::string output_path, std::string condition )
{
  /*~~~~~~~~~~~~~~~~~~*/
  /* 1) Open files    */
  /*~~~~~~~~~~~~~~~~~~*/
  std::stringstream state_trajectory_filename;
  std::stringstream f_trajectory_filename;
  std::stringstream c_trajectory_filename;
  std::stringstream v_trajectory_filename;
  std::stringstream p_trajectory_filename;
  std::stringstream b_trajectory_filename;
  state_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_state_trajectory.csv";
  f_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_f_trajectory.csv";
  c_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_c_trajectory.csv";
  v_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_v_trajectory.csv";
  p_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_p_trajectory.csv";
  b_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_b_trajectory.csv";
  _state_trajectory_file.open(state_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _f_trajectory_file.open(f_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _c_trajectory_file.open(c_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _v_trajectory_file.open(v_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _p_trajectory_file.open(p_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _b_trajectory_file.open(b_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  /*~~~~~~~~~~~~~~~~~~*/
  /* 2) Write headers */
  /*~~~~~~~~~~~~~~~~~~*/
  _state_trajectory_file << "condition;iter;t;dt;mu;doubling_time;density;consistent;mu_diff\n";
  _f_trajectory_file << "condition;iter;t;dt";
  _c_trajectory_file << "condition;iter;t;dt";
  _v_trajectory_file << "condition;iter;t;dt";
  _p_trajectory_file << "condition;iter;t;dt";
  _b_trajectory_file << "condition;iter;t;dt";
  for (int i = 0; i < _nc; i++)
  {
    _c_trajectory_file << ";" << _c_ids[i];
    _b_trajectory_file << ";" << _c_ids[i];
  }
  for (int j = 0; j < _nj; j++)
  {
    _f_trajectory_file << ";" << _reaction_ids[j];
    _v_trajectory_file << ";" << _reaction_ids[j];
    _p_trajectory_file << ";" << _reaction_ids[j];
  }
  _f_trajectory_file << "\n";
  _c_trajectory_file << "\n";
  _v_trajectory_file << "\n";
  _p_trajectory_file << "\n";
  _b_trajectory_file << "\n";
}

/**
 * \brief    Write data into trajectory output files
 * \details  --
 * \param    std::string condition
 * \param    int iter
 * \param    double t
 * \param    double dt
 * \return   \e void
 */
void Model::write_trajectory_output_files( std::string condition, int iter, double t, double dt )
{
  /*------------------------------------*/
  /* 1) Update state file               */
  /*------------------------------------*/
  _state_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt << ";" << _mu << ";" << _doubling_time << ";" << _density << ";" << _consistent << ";" << _mu_diff << "\n";
  _state_trajectory_file.flush();
  /*------------------------------------*/
  /* 2) Update metabolites related file */
  /*------------------------------------*/
  _c_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt;
  _b_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt;
  for (int i = 0; i < _nc; i++)
  {
    _c_trajectory_file << ";" << gsl_vector_get(_c, i);
    _b_trajectory_file << ";" << gsl_vector_get(_b, i);
  }
  _c_trajectory_file << "\n";
  _b_trajectory_file << "\n";
  _c_trajectory_file.flush();
  _b_trajectory_file.flush();
  /*------------------------------------*/
  /* 2) Update reactions related file   */
  /*------------------------------------*/
  _f_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt;
  _v_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt;
  _p_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt;
  for (int j = 0; j < _nj; j++)
  {
    _f_trajectory_file << ";" << gsl_vector_get(_f, j);
    _v_trajectory_file << ";" << gsl_vector_get(_v, j);
    _p_trajectory_file << ";" << gsl_vector_get(_p, j);
  }
  _f_trajectory_file << "\n";
  _v_trajectory_file << "\n";
  _p_trajectory_file << "\n";
  _f_trajectory_file.flush();
  _v_trajectory_file.flush();
  _p_trajectory_file.flush();
  //system("/usr/local/bin/Rscript plot_trajectory.R > /dev/null &");
}

/**
 * \brief    Close trajectory output files
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::close_trajectory_ouput_files( void )
{
  _state_trajectory_file.close();
  _f_trajectory_file.close();
  _c_trajectory_file.close();
  _v_trajectory_file.close();
  _p_trajectory_file.close();
  _b_trajectory_file.close();
}

/**
 * \brief    Open optimum output files
 * \details  Also writes headers
 * \param    std::string output_path
 * \return   \e void
 */
void Model::open_optimum_output_files( std::string output_path, std::string condition )
{
  /*~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Create filenames */
  /*~~~~~~~~~~~~~~~~~~~~~*/
  std::stringstream state_optimum_filename;
  std::stringstream f_optimum_filename;
  std::stringstream c_optimum_filename;
  std::stringstream v_optimum_filename;
  std::stringstream p_optimum_filename;
  std::stringstream b_optimum_filename;
  state_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_state_optimum.csv";
  f_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_f_optimum.csv";
  c_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_c_optimum.csv";
  v_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_v_optimum.csv";
  p_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_p_optimum.csv";
  b_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_b_optimum.csv";
  /*~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Create files     */
  /*~~~~~~~~~~~~~~~~~~~~~*/
  /*** Open files ***/
  _state_optimum_file.open(state_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _f_optimum_file.open(f_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _c_optimum_file.open(c_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _v_optimum_file.open(v_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _p_optimum_file.open(p_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _b_optimum_file.open(b_optimum_filename.str(), std::ios::out | std::ios::trunc);
  /*** Write headers ***/
  _state_optimum_file << "condition;mu;doubling_time;density;consistent;converged;run_time\n";
  _f_optimum_file << "condition";
  _c_optimum_file << "condition";
  _v_optimum_file << "condition";
  _p_optimum_file << "condition";
  _b_optimum_file << "condition";
  for (int i = 0; i < _nc; i++)
  {
    _c_optimum_file << ";" << _c_ids[i];
    _b_optimum_file << ";" << _c_ids[i];
  }
  for (int j = 0; j < _nj; j++)
  {
    _f_optimum_file << ";" << _reaction_ids[j];
    _v_optimum_file << ";" << _reaction_ids[j];
    _p_optimum_file << ";" << _reaction_ids[j];
  }
  _f_optimum_file << "\n";
  _c_optimum_file << "\n";
  _v_optimum_file << "\n";
  _p_optimum_file << "\n";
  _b_optimum_file << "\n";
}

/**
 * \brief    Write data into optimum output files
 * \details  --
 * \param    std::string condition
 * \param    bool converged
 * \param    double runtime
 * \return   \e void
 */
void Model::write_optimum_output_files( std::string condition, bool converged, double runtime )
{
  _state_optimum_file << condition << ";" << _mu << ";" << _doubling_time << ";" << _density << ";" << _consistent << ";" << converged << ";" << runtime << "\n";
  _f_optimum_file << condition;
  _c_optimum_file << condition;
  _v_optimum_file << condition;
  _p_optimum_file << condition;
  _b_optimum_file << condition;
  for (int i = 0; i < _nc; i++)
  {
    _c_optimum_file << ";" << gsl_vector_get(_c, i);
    _b_optimum_file << ";" << gsl_vector_get(_b, i);
  }
  for (int j = 0; j < _nj; j++)
  {
    _f_optimum_file << ";" << gsl_vector_get(_f, j);
    _v_optimum_file << ";" << gsl_vector_get(_v, j);
    _p_optimum_file << ";" << gsl_vector_get(_p, j);
  }
  _f_optimum_file << "\n";
  _c_optimum_file << "\n";
  _v_optimum_file << "\n";
  _p_optimum_file << "\n";
  _b_optimum_file << "\n";
  _f_optimum_file.flush();
  _c_optimum_file.flush();
  _v_optimum_file.flush();
  _p_optimum_file.flush();
  _b_optimum_file.flush();
}

/**
 * \brief    Close optimum output files
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::close_optimum_ouput_files( void )
{
  _state_optimum_file.close();
  _f_optimum_file.close();
  _c_optimum_file.close();
  _v_optimum_file.close();
  _p_optimum_file.close();
  _b_optimum_file.close();
}

/**
 * \brief    Print the optimum to the standard output
 * \details  --
 * \param    std::string condition
 * \param    bool converged
 * \param    double runtime
 * \return   \e void
 */
void Model::print_to_standard_ouput( std::string condition, bool converged, double runtime )
{
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Print condition */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "CONDITION " << condition << std::endl;
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 2) State output    */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "\ncell state" << std::endl;
  std::cout << "condition\tmu\tdoubling_time\tdensity\tconsistent\tconverged\trun_time" << std::endl;
  std::cout << condition << "\t" << _mu << "\t" << _doubling_time << "\t" << _density << "\t" << _consistent << "\t" << converged << "\t" << runtime << "\n";
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 3) f vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "\nf vector" << std::endl;
  std::cout << "condition";
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << _reaction_ids[j];
  }
  std::cout << std::endl << condition;
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << gsl_vector_get(_f, j);
  }
  std::cout << std::endl;
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 4) v vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "\nv vector" << std::endl;
  std::cout << "condition";
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << _reaction_ids[j];
  }
  std::cout << std::endl << condition;
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << gsl_vector_get(_v, j);
  }
  std::cout << std::endl;
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 5) p vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "\np vector" << std::endl;
  std::cout << "condition";
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << _reaction_ids[j];
  }
  std::cout << std::endl << condition;
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << gsl_vector_get(_p, j);
  }
  std::cout << std::endl;
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 6) b vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "\nb vector" << std::endl;
  std::cout << "condition";
  for (int j = 0; j < _nc; j++)
  {
    std::cout << "\t" << _c_ids[j];
  }
  std::cout << std::endl << condition;
  for (int j = 0; j < _nc; j++)
  {
    std::cout << "\t" << gsl_vector_get(_b, j);
  }
  std::cout << std::endl;
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 7) c vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "\nc vector" << std::endl;
  std::cout << "condition";
  for (int j = 0; j < _nc; j++)
  {
    std::cout << "\t" << _c_ids[j];
  }
  std::cout << std::endl << condition;
  for (int j = 0; j < _nc; j++)
  {
    std::cout << "\t" << gsl_vector_get(_c, j);
  }
  std::cout << std::endl;
  std::cout << std::endl;
}

/**
 * \brief    Load metabolite identifiers
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_metabolite_identifiers( void )
{
  assert(_metabolite_ids.size()==0);
  assert(_x_ids.size()==0);
  assert(_c_ids.size()==0);
  assert(_metabolite_indices.size()==0);
  assert(is_file_exist(_model_path+"/"+_model_name+"/M.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/M.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string id;
  getline(file, line);
  int index = 0;
  while(getline(file, line))
  {
    std::stringstream flux(line.c_str());
    getline(flux, id, ';');
    _metabolite_ids.push_back(id);
    _metabolite_indices[id] = index;
    index++;
    if (id.substr(0, 2) == "x_")
    {
      _x_ids.push_back(id);
    }
    else
    {
      _c_ids.push_back(id);
    }
  }
  file.close();
}

/**
 * \brief    Load reaction identifiers
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_reaction_identifiers( void )
{
  assert(_reaction_ids.size()==0);
  assert(_reaction_indices.size()==0);
  assert(is_file_exist(_model_path+"/"+_model_name+"/M.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/M.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string id;
  getline(file, line);
  std::stringstream flux(line.c_str());
  getline(flux, id, ';');
  int index = 0;
  while(getline(flux, id, ';'))
  {
    _reaction_ids.push_back(id);
    _reaction_indices[id] = index;
    index++;
  }
  file.close();
}

/**
 * \brief    Load vector size
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_vector_sizes( void )
{
  _nx = (int)_x_ids.size();
  _nc = (int)_c_ids.size();
  _ni = (int)_metabolite_ids.size();
  _nj = (int)_reaction_ids.size();
}

/**
 * \brief    Load complete and internal M matrices
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_M( void )
{
  assert(_Mx==NULL);
  assert(_M==NULL);
  _Mx = gsl_matrix_alloc(_ni, _nj);
  _M  = gsl_matrix_alloc(_nc, _nj);
  gsl_matrix_set_zero(_Mx);
  gsl_matrix_set_zero(_M);
  int Mx_row = 0;
  int M_row  = 0;
  int col    = 0;
  assert(is_file_exist(_model_path+"/"+_model_name+"/M.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/M.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string id;
  std::string str_value;
  /*** skip header line ***/
  getline(file, line);
  while(getline(file, line))
  {
    std::stringstream flux(line.c_str());
    /*** get metabolite id ***/
    getline(flux, id, ';');
    col = 0;
    /*** fill the matrices for the given line ***/
    while(getline(flux, str_value, ';'))
    {
      double value = stod(str_value);
      gsl_matrix_set(_Mx, Mx_row, col, value);
      if (id.substr(0, 2) != "x_")
      {
        gsl_matrix_set(_M, M_row, col, value);
      }
      col++;
    }
    Mx_row++;
    if (id.substr(0, 2) != "x_")
    {
      M_row++;
    }
  }
  file.close();
}

/**
 * \brief    Load the KM matrix
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_KM( void )
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Load the complete KM matrix             */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  assert(_KM==NULL);
  _KM = gsl_matrix_alloc(_ni, _nj);
  gsl_matrix_set_zero(_KM);
  int row = 0;
  int col = 0;
  assert(is_file_exist(_model_path+"/"+_model_name+"/KM.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/KM.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string id;
  std::string str_value;
  /*** skip header line ***/
  getline(file, line);
  while(getline(file, line))
  {
    std::stringstream flux(line.c_str());
    /*** get metabolite id ***/
    getline(flux, id, ';');
    col = 0;
    /*** fill the matrices for the given line ***/
    while(getline(flux, str_value, ';'))
    {
      double value = stod(str_value);
      gsl_matrix_set(_KM, row, col, value);
      col++;
    }
    row++;
  }
  file.close();
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Create forward and backward KM matrices */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  _KM_f = gsl_matrix_alloc(_ni, _nj);
  _KM_b = gsl_matrix_alloc(_ni, _nj);
  gsl_matrix_set_zero(_KM_f);
  gsl_matrix_set_zero(_KM_b);
  for (int i = 0; i < _ni; i++)
  {
    for (int j = 0; j < _nj; j++)
    {
      if (gsl_matrix_get(_Mx, i, j) < 0.0)
      {
        gsl_matrix_set(_KM_f, i, j, gsl_matrix_get(_KM, i, j));
      }
      else if (gsl_matrix_get(_Mx, i, j) > 0.0)
      {
        gsl_matrix_set(_KM_b, i, j, gsl_matrix_get(_KM, i, j));
      }
    }
  }
}

/**
 * \brief    Load the KI matrix
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_KI( void )
{
  assert(_KI==NULL);
  _KI = gsl_matrix_alloc(_ni, _nj);
  gsl_matrix_set_zero(_KI);
  if (is_file_exist(_model_path+"/"+_model_name+"/KI.csv"))
  {
    int row = 0;
    int col = 0;
    assert(is_file_exist(_model_path+"/"+_model_name+"/KI.csv"));
    std::ifstream file(_model_path+"/"+_model_name+"/KI.csv", std::ios::in);
    assert(file);
    std::string line;
    std::string id;
    std::string str_value;
    /*** skip header line ***/
    getline(file, line);
    while(getline(file, line))
    {
      std::stringstream flux(line.c_str());
      /*** get metabolite id ***/
      getline(flux, id, ';');
      col = 0;
      /*** fill the matrices for the given line ***/
      while(getline(flux, str_value, ';'))
      {
        double value = stod(str_value);
        gsl_matrix_set(_KI, row, col, value);
        col++;
      }
      row++;
    }
    file.close();
  }
}

/**
 * \brief    Load the KA matrix
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_KA( void )
{
  assert(_KA==NULL);
  _KA = gsl_matrix_alloc(_ni, _nj);
  gsl_matrix_set_zero(_KA);
  if (is_file_exist(_model_path+"/"+_model_name+"/KA.csv"))
  {
    int row = 0;
    int col = 0;
    assert(is_file_exist(_model_path+"/"+_model_name+"/KA.csv"));
    std::ifstream file(_model_path+"/"+_model_name+"/KA.csv", std::ios::in);
    assert(file);
    std::string line;
    std::string id;
    std::string str_value;
    /*** skip header line ***/
    getline(file, line);
    while(getline(file, line))
    {
      std::stringstream flux(line.c_str());
      /*** get metabolite id ***/
      getline(flux, id, ';');
      col = 0;
      /*** fill the matrices for the given line ***/
      while(getline(flux, str_value, ';'))
      {
        double value = stod(str_value);
        gsl_matrix_set(_KA, row, col, value);
        col++;
      }
      row++;
    }
    file.close();
  }
}

/**
 * \brief    Load the KR matrix
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_KR( void )
{
  assert(_KR==NULL);
  _KR = gsl_matrix_alloc(_ni, _nj);
  gsl_matrix_set_zero(_KR);
  if (is_file_exist(_model_path+"/"+_model_name+"/KR.csv"))
  {
    int row = 0;
    int col = 0;
    assert(is_file_exist(_model_path+"/"+_model_name+"/KR.csv"));
    std::ifstream file(_model_path+"/"+_model_name+"/KR.csv", std::ios::in);
    assert(file);
    std::string line;
    std::string id;
    std::string str_value;
    /*** skip header line ***/
    getline(file, line);
    while(getline(file, line))
    {
      std::stringstream flux(line.c_str());
      /*** get metabolite id ***/
      getline(flux, id, ';');
      col = 0;
      /*** fill the matrices for the given line ***/
      while(getline(flux, str_value, ';'))
      {
        double value = stod(str_value);
        gsl_matrix_set(_KR, row, col, value);
        col++;
      }
      row++;
    }
    file.close();
  }
}

/**
 * \brief    Load the kcat vectors
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_kcat( void )
{
  assert(_kcat_f==NULL);
  assert(_kcat_b==NULL);
  _kcat_f = gsl_vector_alloc(_nj);
  _kcat_b = gsl_vector_alloc(_nj);
  gsl_vector_set_zero(_kcat_f);
  gsl_vector_set_zero(_kcat_b);
  int col = 0;
  assert(is_file_exist(_model_path+"/"+_model_name+"/kcat.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/kcat.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string id;
  std::string str_value;
  /*** Skip header line ***/
  getline(file, line);
  /*** Load kcat forward values ***/
  getline(file, line);
  std::stringstream flux1(line.c_str());
  getline(flux1, id, ';');
  col = 0;
  while(getline(flux1, str_value, ';'))
  {
    double value = stod(str_value);
    gsl_vector_set(_kcat_f, col, value);
    col++;
  }
  /*** Load kcat backward values ***/
  getline(file, line);
  std::stringstream flux2(line.c_str());
  getline(flux2, id, ';');
  col = 0;
  while(getline(flux2, str_value, ';'))
  {
    double value = stod(str_value);
    gsl_vector_set(_kcat_b, col, value);
    col++;
  }
  /*** Close the file ***/
  file.close();
}

/**
 * \brief    Load conditions
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_conditions( void )
{
  assert(_condition_ids.size()==0);
  assert(_condition_params.size()==0);
  assert(_condition_indices.size()==0);
  assert(_conditions==NULL);
  assert(is_file_exist(_model_path+"/"+_model_name+"/conditions.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/conditions.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string id;
  std::string str_value;
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Read condition identifiers */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  getline(file, line);
  std::stringstream flux1(line.c_str());
  getline(flux1, id, ';');
  int index = 0;
  while(getline(flux1, id, ';'))
  {
    _condition_ids.push_back(id);
    _condition_indices[id] = index;
    index++;
  }
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Read condition parameters  */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  while(getline(file, line))
  {
    std::stringstream flux(line.c_str());
    getline(flux, id, ';');
    _condition_params.push_back(id);
  }
  file.close();
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Load the conditions        */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  _conditions = gsl_matrix_alloc(_condition_params.size(), _condition_ids.size());
  gsl_matrix_set_zero(_conditions);
  int row = 0;
  file.open(_model_path+"/"+_model_name+"/conditions.csv", std::ios::in);
  assert(file);
  getline(file, line);
  while (getline(file, line))
  {
    std::stringstream flux(line.c_str());
    getline(flux, id, ';');
    int col = 0;
    while(getline(flux, str_value, ';'))
    {
      double value = stod(str_value);
      gsl_matrix_set(_conditions, row, col, value);
      col++;
    }
    row++;
  }
  file.close();
}

/**
 * \brief    Load constant reactions
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_constant_reactions( void )
{
  _constant_reactions.clear();
  assert(is_file_exist(_model_path+"/"+_model_name+"/constant_reactions.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/constant_reactions.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string reaction_id;
  std::string str_value;
  getline(file, line);
  while(getline(file, line))
  {
    std::stringstream flux(line.c_str());
    getline(flux, reaction_id, ';');
    getline(flux, str_value, ';');
    assert(_constant_reactions.find(reaction_id) == _constant_reactions.end());
    _constant_reactions[reaction_id] = stod(str_value);
  }
  file.close();
}

/**
 * \brief    Load f0
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_f0( void )
{
  assert(_f0==NULL);
  _f0 = gsl_vector_alloc(_nj);
  gsl_vector_set_zero(_f0);
  assert(is_file_exist(_model_path+"/"+_model_name+"/f0.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/f0.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string reaction_id;
  std::string str_value;
  getline(file, line);
  while(getline(file, line))
  {
    std::stringstream flux(line.c_str());
    getline(flux, reaction_id, ';');
    getline(flux, str_value, ';');
    assert(_reaction_indices.find(reaction_id) != _reaction_indices.end());
    double value = stod(str_value);
    gsl_vector_set(_f0, _reaction_indices[reaction_id], value);
  }
  file.close();
}

/**
 * \brief    Re-load f0
 * \details  Re-load the vector from the last trajectory point
 * \param    void
 * \return   \e void
 */
void Model::reload_f0( void )
{
  assert(_f0==NULL);
  _f0 = gsl_vector_alloc(_nj);
  gsl_vector_set_zero(_f0);
  assert(is_file_exist(_model_path+"/"+_model_name+"/f0.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/f0.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string reaction_id;
  std::string str_value;
  getline(file, line);
  while(getline(file, line))
  {
    std::stringstream flux(line.c_str());
    getline(flux, reaction_id, ';');
    getline(flux, str_value, ';');
    assert(_reaction_indices.find(reaction_id) != _reaction_indices.end());
    double value = stod(str_value);
    gsl_vector_set(_f0, _reaction_indices[reaction_id], value);
  }
  file.close();
}

/**
 * \brief    Initialize all mathematical variables
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::initialize_variables( void )
{
  initialize_static_variables();
  initialize_dynamic_variables();
}

/**
 * \brief    Initialize static variables
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::initialize_static_variables( void )
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Initialize constants      */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  assert(_sM==NULL);
  _sM = gsl_vector_alloc(_nj);
  gsl_vector_set_zero(_sM);
  for(int i = 0; i < _nc; i++)
  {
    for(int j = 0; j < _nj; j++)
    {
      double current = gsl_vector_get(_sM, j);
      gsl_vector_set(_sM, j, current+gsl_matrix_get(_M, i, j));
    }
  }
  _r = _nj-1;
  _a = _nc-1;
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Initialize reaction types */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  _type              = new rtype[_nj];
  gsl_vector* KI_vec = gsl_vector_alloc(_ni);
  gsl_vector* KA_vec = gsl_vector_alloc(_ni);
  gsl_vector* KR_vec = gsl_vector_alloc(_ni);
  for(int j = 0; j < _nj; j++)
  {
    gsl_matrix_get_col(KI_vec, _KI, j);
    gsl_matrix_get_col(KA_vec, _KA, j);
    gsl_matrix_get_col(KR_vec, _KR, j);
    bool kcat_b_zero = fabs(gsl_vector_get(_kcat_b, j)) < _tol;
    bool KI_sum_zero = fabs(gsl_blas_dasum(KI_vec)) < _tol;
    bool KA_sum_zero = fabs(gsl_blas_dasum(KA_vec)) < _tol;
    bool KR_sum_zero = fabs(gsl_blas_dasum(KR_vec)) < _tol;
    if (kcat_b_zero && KI_sum_zero && KA_sum_zero && KR_sum_zero)
    {
      _type[j] = IMM;
    }
    else if (kcat_b_zero && !KI_sum_zero && KA_sum_zero && KR_sum_zero)
    {
      _type[j] = IMMI;
    }
    else if (kcat_b_zero && KI_sum_zero && !KA_sum_zero && KR_sum_zero)
    {
      _type[j] = IMMA;
    }
    else if (kcat_b_zero && !KI_sum_zero && !KA_sum_zero && KR_sum_zero)
    {
      _type[j] = IMMIA;
    }
    else if (!kcat_b_zero && KR_sum_zero)
    {
      assert(KI_sum_zero);
      assert(KA_sum_zero);
      _type[j] = RMM;
    }
    else if (kcat_b_zero && !KR_sum_zero)
    {
      assert(KI_sum_zero);
      assert(KA_sum_zero);
      _type[j] = IMMR;
    }
  }
  gsl_vector_free(KI_vec);
  gsl_vector_free(KA_vec);
  gsl_vector_free(KR_vec);
  KI_vec = NULL;
  KA_vec = NULL;
  KR_vec = NULL;
}

/**
 * \brief    Initialize dynamic variables
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::initialize_dynamic_variables( void )
{
  /*** Assertions ***/
  assert(_x==NULL);
  assert(_f==NULL);
  assert(_f_trunc==NULL);
  assert(_c==NULL);
  assert(_xc==NULL);
  assert(_tau_j==NULL);
  assert(_v==NULL);
  assert(_p==NULL);
  assert(_b==NULL);
  assert(_ditau_j==NULL);
  assert(_dmu_f==NULL);
  assert(_GCC_f==NULL);
  assert(_dmu_f_term2==NULL);
  assert(_dmu_f_term3==NULL);
  assert(_dmu_f_term4==NULL);
  assert(_dmu_f_term5==NULL);
  /*** Allocate memory ***/
  _x           = gsl_vector_alloc(_nx);
  _f           = gsl_vector_alloc(_nj);
  _f_trunc     = gsl_vector_alloc(_nj-1);
  _c           = gsl_vector_alloc(_nc);
  _xc          = gsl_vector_alloc(_ni);
  _tau_j       = gsl_vector_alloc(_nj);
  _v           = gsl_vector_alloc(_nj);
  _p           = gsl_vector_alloc(_nj);
  _b           = gsl_vector_alloc(_nc);
  _ditau_j     = gsl_matrix_alloc(_nj, _nc);
  _dmu_f       = gsl_vector_alloc(_nj);
  _GCC_f       = gsl_vector_alloc(_nj);
  _dmu_f_term2 = gsl_vector_alloc(_nj);
  _dmu_f_term3 = gsl_matrix_alloc(_nj, _nj);
  _dmu_f_term4 = gsl_vector_alloc(_nj);
  _dmu_f_term5 = gsl_vector_alloc(_nj);
  /*** Initialize all variables to zero ***/
  gsl_vector_set_zero(_x);
  gsl_vector_set_zero(_f);
  gsl_vector_set_zero(_f_trunc);
  gsl_vector_set_zero(_c);
  gsl_vector_set_zero(_xc);
  gsl_vector_set_zero(_tau_j);
  gsl_vector_set_zero(_v);
  gsl_vector_set_zero(_p);
  gsl_vector_set_zero(_b);
  gsl_matrix_set_zero(_ditau_j);
  gsl_vector_set_zero(_dmu_f);
  gsl_vector_set_zero(_GCC_f);
  gsl_vector_set_zero(_dmu_f_term2);
  gsl_matrix_set_zero(_dmu_f_term3);
  gsl_vector_set_zero(_dmu_f_term4);
  gsl_vector_set_zero(_dmu_f_term5);
  /*** Initialize vector views ***/
  _x_view = gsl_vector_subvector(_xc, 0, _nx);
  _c_view = gsl_vector_subvector(_xc, _nx, _nc);
}

/**
 * \brief    Return a Gaussian kernel term
 * \details  --
 * \param    double x
 * \param    double mu
 * \return   \e double
 */
double Model::gaussian_term( double x, double mu )
{
  double sigma = REGULATION_SIGMA*x;
  double val = (x-mu)/(sigma*sigma);
  return(val);
}

/**
 * \brief    Return a Gaussian kernel value
 * \details  --
 * \param    double x
 * \param    double mu
 * \return   \e double
 */
double Model::gaussian_kernel( double x, double mu )
{
  double sigma = REGULATION_SIGMA*x;
  double val = (x-mu)/sigma;
  double res = -0.5*val*val;
  if (res < -700)
  {
    return(0.0);
  }
  else
  {
    return(gsl_sf_exp(res));
  }
}

/**
 * \brief    Calculate model state
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::calculate( void )
{
  calculate_first_order_terms();
  calculate_second_order_terms();
  check_model_consistency();
}

/**
 * \brief    Compute internal concentrations
 * \details  Formula: rho*M*f
 * \param    void
 * \return   \e void
 */
void Model::compute_c( void )
{
  gsl_blas_dgemv(CblasNoTrans, _current_rho, _M, _f, 0.0, _c);
  if (_adjust_concentrations)
  {
    for (int i = 0; i < _nc; i++)
    {
      if (gsl_vector_get(_c, i) < _tol)
      {
        gsl_vector_set(_c, i, _tol);
      }
    }
  }
}

/**
 * \brief    Compute the concentration vector
 * \details  Simple concatenation
 * \param    void
 * \return   \e void
 */
void Model::compute_xc( void )
{
  gsl_vector_memcpy(&_x_view.vector, _x);
  gsl_vector_memcpy(&_c_view.vector, _c);
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics
 * \details  Formula: tau_j = prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMM( int j )
{
  double term1 = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    term1 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  double term2 = gsl_vector_get(_kcat_f, j);
  gsl_vector_set(_tau_j, j, term1/term2);
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + inhibition (only one inhibitor per reaction)
 * \details  Formula: tau_j = prod(1+xc*1/KI[,j])*prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMMi( int j )
{
  double term1 = 1.0;
  double term2 = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    double rKI = (gsl_matrix_get(_KI, i, j) > 0.0 ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    term1     *= 1.0+gsl_vector_get(_xc, i)*rKI;
    term2     *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  double term3 = gsl_vector_get(_kcat_f, j);
  gsl_vector_set(_tau_j, j, term1*term2/term3);
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + activation (only one activator per reaction)
 * \details  Formula: tau_j = prod(1+KA[,j]/xc)*prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMMa( int j )
{
  double term1 = 1.0;
  double term2 = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    term1 *= 1.0+gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
    term2 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  double term3 = gsl_vector_get(_kcat_f, j);
  gsl_vector_set(_tau_j, j, term1*term2/term3);
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + inhibition + activation
 * \details  Formula: tau_j = prod(1+xc*1/KI[,j])*prod(1+KA[,j]/xc)*prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMMia( int j )
{
  double term1 = 1.0;
  double term2 = 1.0;
  double term3 = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    double rKI = (gsl_matrix_get(_KI, i, j) > 0.0 ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    term1     *= 1.0+gsl_vector_get(_xc, i)*rKI;
    term2     *= 1.0+gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
    term3     *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  double term4 = gsl_vector_get(_kcat_f, j);
  gsl_vector_set(_tau_j, j, term1*term2*term3/term4);
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + regulation
 * \details  Formula: tau_j = prod((Km_f[,j]+xc)/(xc*f(xc))) * 1/kcat_f[j] with f(xc) = exp( -0.5*((xc-JR[:j])/10.0)^2 )
 * \param    int j
 * \return   \e void
 */
void Model::iMMr( int j )
{
  double term1 = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    double x       = gsl_vector_get(_xc, i);
    double KM      = gsl_matrix_get(_KM_f, i, j);
    double KR      = gsl_matrix_get(_KR, i, j) > _tol ? gsl_matrix_get(_KR, i, j) : x;
    double gkernel = gaussian_kernel(x, KR);
    term1         *= (KM+x)/(x*gkernel);
  }
  double term2 = gsl_vector_get(_kcat_f, j);
  gsl_vector_set(_tau_j, j, term1/term2);
}

/**
 * \brief    Reversible Michaelis-Menten kinetics
 * \details  Formula: tau_j = 1/[ kcat_f[j]/prod(1+Km_f[,j]/xc) - kcat_b[j]/prod(1+Km_b[,j]/xc) ]
 * \param    int j
 * \return   \e void
 */
void Model::rMM( int j )
{
  double term1 = gsl_vector_get(_kcat_f, j);
  double term2 = 1.0;
  double term3 = gsl_vector_get(_kcat_b, j);
  double term4 = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    term2 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    term4 *= 1.0+gsl_matrix_get(_KM_b, i, j)/gsl_vector_get(_xc, i);
  }
  gsl_vector_set(_tau_j, j, 1.0/(term1/term2-term3/term4));
}

/**
 * \brief    Compute tau_j
 * \details  --
 * \param    int j
 * \return   \e void
 */
void Model::compute_tau( int j )
{
  switch(_type[j])
  {
    case IMM:
      iMM(j);
      break;
    case IMMI:
      iMMi(j);
      break;
    case IMMA:
      iMMa(j);
      break;
    case IMMIA:
      iMMia(j);
      break;
    case IMMR:
      iMMr(j);
    case RMM:
      rMM(j);
      break;
  }
}

/**
 * \brief    Derivative of iMM with respect to metabolite concentrations
 * \details  Formula: --
 * \param    int j
 * \return   \e void
 */
void Model::diMM( int j )
{
  double constant1 = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double term1 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term2 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term2 *= 1.0+gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
      gsl_matrix_set(_ditau_j, j, i, -term1*term2/constant1);
    }
  }
}

/**
 * \brief    Derivative of iMMi with respect to metabolite concentrations
 * \details  Formula: --
 * \param    int j
 * \return   \e void
 */
void Model::diMMi( int j )
{
  double constant1 = 1.0;
  double constant2 = 1.0;
  double constant3 = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    double rKI = (gsl_matrix_get(_KI, i, j) > 1e-10 ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    constant1 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    constant2 *= 1.0+gsl_vector_get(_xc, i)*rKI;
  }
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double rKI   = (gsl_matrix_get(_KI, y, j) > 1e-10 ? 1.0/gsl_matrix_get(_KI, y, j) : 0.0);
    double term1 = rKI*constant1;
    double term2 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term3 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term3 *= 1.0+gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
    }
    gsl_matrix_set(_ditau_j, j, i, (term1-constant2*term2*term3)/constant3);
  }
}

/**
 * \brief    Derivative of iMMa with respect to metabolite concentrations
 * \details  Formula: --
 * \param    int j
 * \return   \e void
 */
void Model::diMMa( int j )
{
  double constant1 = 1.0;
  double constant2 = 1.0;
  double constant3 = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    constant1 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    constant2 *= 1.0+gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
  }
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double term1 = gsl_matrix_get(_KA, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term2 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term3 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term3 *= 1.0+gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
    }
    gsl_matrix_set(_ditau_j, j, i, -(constant1*term1+constant2*term2*term3)/constant3);
  }
}

/**
 * \brief    Derivative of iMMia with respect to metabolite concentrations
 * \details  Formula: --
 * \param    int j
 * \return   \e void
 */
void Model::diMMia( int j )
{
  double constant1 = 1.0;
  double constant2 = 1.0;
  double constant3 = 1.0;
  double constant4 = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    double rKI = (gsl_matrix_get(_KI, i, j) > 0.0 ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    constant1 *= 1.0+gsl_vector_get(_xc, i)*rKI;
    constant2 *= 1.0+gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
    constant3 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double rKI   = (gsl_matrix_get(_KI, y, j) > 0.0 ? 1.0/gsl_matrix_get(_KI, y, j) : 0.0);
    double term1 = rKI;
    double term2 = -gsl_matrix_get(_KA, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term3 = -gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term4 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term4 *= 1.0+gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
    }
    double term5 = (term1*constant2*constant3)+(term2*constant1*constant3)+(term3*term4*constant1*constant2);
    gsl_matrix_set(_ditau_j, j, i, term5/constant4);
  }
}

/**
 * \brief    Derivative of iMMr with respect to metabolite concentrations
 * \details  Formula: --
 * \param    int j
 * \return   \e void
 */
void Model::diMMr( int j )
{
  double constant1 = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _nc; i++)
  {
    int    y       = i+_nx;
    double x       = gsl_vector_get(_c, i);
    double KM      = gsl_matrix_get(_KM_f, y, j);
    double KR      = gsl_matrix_get(_KR, y, j) > _tol ? gsl_matrix_get(_KR, y, j) : x;
    double gkernel = gaussian_kernel(x, KR);
    double gterm   = gaussian_term(x, KR);
    double term1   = gkernel*(-KM/(x*x)+(KM+x)/x*gterm);
    double term2   = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        double x       = gsl_vector_get(_xc, index);
        double KR      = gsl_matrix_get(_KR, index, j) > _tol ? gsl_matrix_get(_KR, index, j) : x;
        double gkernel = gaussian_kernel(x, KR);
        term2         *= (gsl_matrix_get(_KM_f, index, j)+x)/(x*gkernel);
      }
    }
    gsl_matrix_set(_ditau_j, j, i, term1*term2/constant1);
  }
}

/**
 * \brief    Derivative of rMM with respect to metabolite concentrations
 * \details  Formula: --
 * \param    int j
 * \return   \e void
 */
void Model::drMM( int j )
{
  double constant1 = gsl_vector_get(_kcat_f, j);
  double constant2 = gsl_vector_get(_kcat_b, j);
  double constant3 = 1.0;
  double constant4 = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    constant3 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    constant4 *= 1.0+gsl_matrix_get(_KM_b, i, j)/gsl_vector_get(_xc, i);
  }
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double term1 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i)+gsl_matrix_get(_KM_f, y, j), 2);
    double term2 = gsl_matrix_get(_KM_b, y, j)/gsl_pow_int(gsl_vector_get(_c, i)+gsl_matrix_get(_KM_b, y, j), 2);
    double term3 = 1.0;
    double term4 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term1 *= 1.0+gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
        term3 *= 1.0+gsl_matrix_get(_KM_b, index, j)/gsl_vector_get(_xc, index);
      }
    }
    double term5 = term1*constant1/term3-term2*constant2/term4;
    double term6 = constant1/constant3-constant2/constant4;
    gsl_matrix_set(_ditau_j, j, i, -term5/gsl_pow_int(term6, 2));
  }
}

/**
 * \brief    Compute ditau_j
 * \details  --
 * \param    int j
 * \return   \e void
 */
void Model::compute_dtau( int j )
{
  switch(_type[j])
  {
    case IMM:
      diMM(j);
      break;
    case IMMI:
      diMMi(j);
      break;
    case IMMA:
      diMMa(j);
      break;
    case IMMIA:
      diMMia(j);
      break;
    case IMMR:
      diMMr(j);
      break;
    case RMM:
      drMM(j);
      break;
  }
}

/**
 * \brief    Compute the growth rate mu
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::compute_mu( void )
{
  double term1 = 0.0;
  gsl_blas_ddot(_tau_j, _f, &term1);
  _mu            = gsl_matrix_get(_M, _a, _r)*gsl_vector_get(_f, _r)/term1;
  _doubling_time = std::log(2)/std::log(1.0 + _mu);
}

/**
 * \brief    Compute the mass flux vector v
 * \details  v = mu*rho*f
 * \param    void
 * \return   \e void
 */
void Model::compute_v( void )
{
  gsl_vector_memcpy(_v, _f);
  gsl_vector_scale(_v, _mu*_current_rho);
}

/**
 * \brief    Compute protein mass concentrations p
 * \details  p = tau_j*v
 * \param    void
 * \return   \e void
 */
void Model::compute_p( void )
{
  gsl_vector_memcpy(_p, _v);
  gsl_vector_mul(_p, _tau_j);
}

/**
 * \brief    Compute mass fractions b
 * \details  b = M*f
 * \param    void
 * \return   \e void
 */
void Model::compute_b( void )
{
  gsl_blas_dgemv(CblasNoTrans, 1.0, _M, _f, 0.0, _b);
}

/**
 * \brief    Compute cell density (should be always 1)
 * \details  density = sM*f
 * \param    void
 * \return   \e void
 */
void Model::compute_density( void )
{
  gsl_blas_ddot(_sM, _f, &_density);
}

/**
 * \brief    Compute local mu gradient with respect to f
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::compute_dmu_f( void )
{
  /*--------*/
  _dmu_f_term1 = gsl_pow_int(_mu, 2)/gsl_vector_get(_b, _a);
  /*--------*/
  gsl_matrix_get_row(_dmu_f_term2, _M, _a);
  gsl_vector_scale(_dmu_f_term2, 1.0/_mu);
  /*--------*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, _current_rho, _ditau_j, _M, 0.0, _dmu_f_term3);
  gsl_blas_dgemv(CblasTrans, 1.0, _dmu_f_term3, _f, 0.0, _dmu_f_term4);
  /*--------*/
  gsl_vector_memcpy(_dmu_f_term5, _tau_j);
  /*--------*/
  gsl_vector_memcpy(_dmu_f, _dmu_f_term2);
  gsl_vector_sub(_dmu_f, _dmu_f_term4);
  gsl_vector_sub(_dmu_f, _dmu_f_term5);
  gsl_vector_scale(_dmu_f, _dmu_f_term1);
}

/**
 * \brief    Compute local growth control coefficients with respect to f
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::compute_GCC_f( void )
{
  // self.dmu_f-self.dmu_f[0]*(self.sM/self.sM[0])
  gsl_vector_memcpy(_GCC_f, _sM);
  gsl_vector_scale(_GCC_f, gsl_vector_get(_dmu_f, 0)/gsl_vector_get(_sM, 0));
  gsl_vector_scale(_GCC_f, -1.0);
  gsl_vector_add(_GCC_f, _dmu_f);
}

/**
 * \brief    Calculate all first order variables from the f vector
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::calculate_first_order_terms( void )
{
  compute_c();
  compute_xc();
  for (int j = 0; j < _nj; j++)
  {
    compute_tau(j);
  }
  compute_mu();
  compute_v();
  compute_p();
  compute_b();
  compute_density();
}

/**
 * \brief    Calculate all secon order variables from the f vector
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::calculate_second_order_terms( void )
{
  for (int j = 0; j < _nj; j++)
  {
    compute_dtau(j);
  }
  compute_dmu_f();
  compute_GCC_f();
}

/**
 * \brief    Check model consistency
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::check_model_consistency( void )
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Test density constraint                 */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  bool test1 = true;
  if (fabs(_density-1.0) >= _tol)
  {
    test1 = false;
  }
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Test negative concentrations constraint */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  bool test2 = true;
  if (gsl_vector_min(_c) < -_tol)
  {
    test2 = false;
  }
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Test negative proteins constraint       */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  bool test3 = true;
  if (gsl_vector_min(_p) < -_tol)
  {
    test3 = false;
  }
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) Print error message if inconsistent     */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  _consistent = true;
  if (!(test1 && test2 && test3))
  {
    _consistent = false;
  }
}

/**
 * \brief    Detect reactions to block
 * \details  Reactions tending to zero or enforcing reversibility
 * \param    void
 * \return   \e void
 */
void Model::block_reactions( void )
{
  for (int j = 0; j < _nj-1; j++)
  {
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* 1) Reaction is irreversible and positive    */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if (_type[j+1] != RMM && gsl_vector_get(_f_trunc, j) <= _tol)
    {
      gsl_vector_set(_f_trunc, j, _tol);
      if (gsl_vector_get(_GCC_f, j+1) < 0.0)
      {
        gsl_vector_set(_GCC_f, j+1, 0.0);
      }
    }
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* 2) Reaction is irreversible and negative    */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if (_type[j+1] == RMM && gsl_vector_get(_f_trunc, j) >= -_tol)
    {
      gsl_vector_set(_f_trunc, j, -_tol);
      if (gsl_vector_get(_GCC_f, j+1) > 0.0)
      {
        gsl_vector_set(_GCC_f, j+1, 0.0);
      }
    }
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* 3) Reaction is reversible and tends to zero */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    else if (_type[j+1] == RMM && fabs(gsl_vector_get(_f_trunc, j)) <= _tol)
    {
      gsl_vector_set(_GCC_f, j+1, 0.0);
      if (gsl_vector_get(_f_trunc, j) >= 0.0)
      {
        gsl_vector_set(_f_trunc, j, _tol);
      }
      else if (gsl_vector_get(_f_trunc, j) < 0.0)
      {
        gsl_vector_set(_f_trunc, j, -_tol);
      }
    }
  }
  for (auto item : _constant_reactions)
  {
    int    j   = _reaction_indices[item.first];
    double val = item.second;
    assert(j > 0);
    gsl_vector_set(_f_trunc, j-1, val);
    gsl_vector_set(_GCC_f, j, 0.0);
  }
}

