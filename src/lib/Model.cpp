/**
 * \file      Model.cpp
 * \author    Charles Rocabert
 * \date      22-07-2024
 * \copyright gbacpp. Copyright © 2024-2025 Charles Rocabert.
 * \license   GNU General Public License v3 (GPLv3)
 * \brief     Model class definition
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
  
  _tol    = 1e-10;
  _mu_tol = 1e-10;
  
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
  _K          = NULL;
  _KM_f       = NULL;
  _KM_b       = NULL;
  _KI         = NULL;
  _KA         = NULL;
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
  
  _condition = "";
  _rho       = 0.0;
  _x         = NULL;
  
  /*----------------------------------------------- GBA first order variables */
  
  _q0                    = NULL;
  _q                     = NULL;
  _q_trunc               = NULL;
  _c                     = NULL;
  _xc                    = NULL;
  _tau_j                 = NULL;
  _v                     = NULL;
  _p                     = NULL;
  _b                     = NULL;
  _density               = 0.0;
  _mu                    = 0.0;
  _consistent            = false;
  _adjust_concentrations = false;
  
  /*----------------------------------------------- GBA second order variables */
  
  _ditau_j = NULL;
  _dmu_dq  = NULL;
  _Gamma   = NULL;
  
  /*----------------------------------------------- Variables for calculation and optimization */
  
  _dmu_dq_term1  = 0.0;
  _dmu_dq_term2  = NULL;
  _dmu_dq_term3  = NULL;
  _dmu_dq_term4  = NULL;
  _dmu_dq_term5  = NULL;
  _stable_count  = 0;
  _mu_diff       = 0.0;
  _mu_rel_diff   = 0.0;
  
  /*----------------------------------------------- Solutions */
  
  _nb_random_solutions = 0;
  _random_solutions.clear();
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
  gsl_matrix_free(_K);
  gsl_matrix_free(_KM_f);
  gsl_matrix_free(_KM_b);
  gsl_matrix_free(_KI);
  gsl_matrix_free(_KA);
  gsl_vector_free(_kcat_f);
  gsl_vector_free(_kcat_b);
  delete[] _type;
  gsl_matrix_free(_conditions);
  _Mx         = NULL;
  _M          = NULL;
  _K          = NULL;
  _KM_f       = NULL;
  _KM_b       = NULL;
  _KI         = NULL;
  _KA         = NULL;
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
  
  gsl_vector_free(_q0);
  gsl_vector_free(_q);
  gsl_vector_free(_q_trunc);
  gsl_vector_free(_c);
  gsl_vector_free(_xc);
  gsl_vector_free(_tau_j);
  gsl_vector_free(_v);
  gsl_vector_free(_p);
  gsl_vector_free(_b);
  _q0         = NULL;
  _q          = NULL;
  _q_trunc    = NULL;
  _c          = NULL;
  _xc         = NULL;
  _tau_j      = NULL;
  _v          = NULL;
  _p          = NULL;
  _b          = NULL;
  
  /*----------------------------------------------- GBA second order variables */
  
  gsl_matrix_free(_ditau_j);
  gsl_vector_free(_dmu_dq);
  gsl_vector_free(_Gamma);
  _ditau_j = NULL;
  _dmu_dq  = NULL;
  _Gamma   = NULL;
  
  /*----------------------------------------------- Variables for calculation optimization */
  
  gsl_vector_free(_dmu_dq_term2);
  gsl_matrix_free(_dmu_dq_term3);
  gsl_vector_free(_dmu_dq_term4);
  gsl_vector_free(_dmu_dq_term5);
  _dmu_dq_term2 = NULL;
  _dmu_dq_term3 = NULL;
  _dmu_dq_term4 = NULL;
  _dmu_dq_term5 = NULL;
  
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
  load_K();
  load_KI();
  load_KA();
  load_kcat();
  load_conditions();
  load_constant_reactions();
  load_q0();
  initialize_variables();
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
 * \param    bool write_optimum
 * \param    bool write_trajectory
 * \param    std::string output_path
 * \param    int stable_count
 * \param    int max_iter
 * \param    bool hessian
 * \param    bool reload
 * \param    bool restart
 * \param    bool verbose
 * \param    bool extra_verbose
 * \return   \e bool
 */
void Model::compute_optimum( std::string condition, bool print_optimum, bool write_optimum, bool write_trajectory, std::string output_path, int stable_count, int max_iter, bool hessian, bool reload, bool restart, bool verbose, bool extra_verbose )
{
  std::clock_t begin = clock();
  bool converged     = compute_gradient_ascent(condition, write_trajectory, output_path, stable_count, max_iter, hessian, reload, restart, verbose, extra_verbose);
  std::clock_t end   = clock();
  double runtime     = double(end-begin)/CLOCKS_PER_SEC;
  if (write_optimum)
  {
    open_optimum_output_files(output_path, condition);
    write_optimum_output_files(condition, converged, runtime);
    close_optimum_ouput_files();
  }
  if (print_optimum)
  {
    print_to_standard_ouput(condition, converged, runtime);
  }
  if ((verbose || extra_verbose) && converged)
  {
    std::cout << " > Condition " << condition << ": convergence reached (mu=" << _mu << ", runtime=" << runtime << ")" << std::endl;
  }
  else if ((verbose || extra_verbose) && !converged)
  {
    std::cout << " > Condition " << condition << ": convergence not reached after " << max_iter << " iterations (mu=" << _mu << ", runtime=" << runtime << ")" << std::endl;
  }
}

/**
 * \brief    Compute the optimum for all the conditions
 * \details  --
 * \param    bool print_optimum
 * \param    bool write_optimum
 * \param    bool write_trajectory
 * \param    std::string output_path
 * \param    int stable_count
 * \param    int max_iter
 * \param    bool hessian
 * \param    bool reload
 * \param    bool restart
 * \param    bool use_previous_sol
 * \param    bool verbose
 * \param    bool extra_verbose
 * \return   \e bool
 */
void Model::compute_optimum_by_condition( bool print_optimum, bool write_optimum, bool write_trajectory, std::string output_path, int stable_count, int max_iter, bool hessian, bool reload, bool restart, bool use_previous_sol, bool verbose, bool extra_verbose )
{
  if (write_optimum)
  {
    open_optimum_output_files(output_path, "all");
  }
  bool reload_local  = false;
  bool restart_local = false;
  for (int i = 0; i < (int)_condition_ids.size(); i++)
  {
    if (i > 0 && use_previous_sol)
    {
      reload_local  = true;
      restart_local = true;
      std::stringstream source;
      std::stringstream target;
      std::stringstream cmdline;
      if (i < (int)_condition_ids.size())
      {
        source << output_path << "/" << _model_name << "_" << _condition_ids[(i-1)] << "_q.bin";
        target << output_path << "/" << _model_name << "_" <<  _condition_ids[i] << "_q.bin";
        cmdline << "cp " << source.str() << " " << target.str();
        std::system(cmdline.str().c_str());
      }
    }
    std::clock_t begin     = clock();
    std::string  condition = _condition_ids[i];
    bool         converged = compute_gradient_ascent(condition, write_trajectory, output_path, stable_count, max_iter, hessian, reload_local, restart_local, verbose, extra_verbose);
    std::clock_t end       = clock();
    double       runtime   = double(end-begin)/CLOCKS_PER_SEC;
    if (write_optimum)
    {
      write_optimum_output_files(condition, converged, runtime);
    }
    if (print_optimum)
    {
      print_to_standard_ouput(condition, converged, runtime);
    }
    if ((verbose || extra_verbose) && converged)
    {
      std::cout << " > Condition " << condition << ": convergence reached (mu=" << _mu << ", runtime=" << runtime << ")" << std::endl;
    }
    else if ((verbose || extra_verbose) && !converged)
    {
      std::cout << " > Condition " << condition << ": convergence not reached after " << max_iter << " iterations (mu=" << _mu << ", runtime=" << runtime << ")" << std::endl;
    }
  }
  if (write_optimum)
  {
    close_optimum_ouput_files();
  }
}

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
 * \param    int max_iter
 * \param    bool hessian
 * \param    bool reload
 * \param    bool restart
 * \param    bool verbose
 * \param    bool extra_verbose
 * \return   \e bool
 */
bool Model::compute_gradient_ascent( std::string condition, bool write_trajectory, std::string output_path, int stable_count, int max_iter, bool hessian, bool reload, bool restart, bool verbose, bool extra_verbose )
{
  auto it = std::find(_condition_ids.begin(), _condition_ids.end(), condition);
  if (it==_condition_ids.end())
  {
    throw std::invalid_argument("> Error: Unknown condition");
  }
  if (!is_path_exist(output_path))
  {
    throw std::invalid_argument("> Error: Path '"+output_path+"' does not exist");
  }
  if (stable_count < 0)
  {
    throw std::invalid_argument("> Error: The stable count parameter must be positive or null");
  }
  if (max_iter <= 0)
  {
    throw std::invalid_argument("> Error: The maximum number of iterations must be positive");
  }
  if (write_trajectory)
  {
    bool append = reload && !restart;
    open_trajectory_output_files(output_path, condition, append);
  }
  double previous_mu   = 0.0;
  double t             = 0.0;
  double dt            = 0.01;
  int    dt_counter    = 0;
  int    nb_iterations = 0;
  _stable_count        = 0;
  _mu_diff             = 0.0;
  _mu_rel_diff         = 0.0;
  if (reload)
  {
    reload_q0(nb_iterations, t, dt, output_path, condition, restart);
  }
  _adjust_concentrations = false;
  set_condition(condition);
  initialize_q();
  calculate();
  if (!_consistent)
  {
    throw std::runtime_error("> Error: The initial solution q0 is not consistent");
  }
  if (extra_verbose)
  {
    std::cout << " > Initial growth rate = " << _mu << std::endl;
  }
  gsl_vector* previous_q_trunc         = gsl_vector_alloc(_nj-1);
  gsl_vector* hessian_previous_q_trunc = gsl_vector_alloc(_nj-1);
  gsl_vector* previous_Gamma_trunc     = gsl_vector_alloc(_nj-1);
  gsl_vector* scaled_Gammadt_trunc     = gsl_vector_alloc(_nj-1);
  gsl_vector_view Gamma_trunc          = gsl_vector_subvector(_Gamma, 1, _nj-1);
  gsl_vector_memcpy(previous_q_trunc, _q_trunc);
  gsl_vector_memcpy(hessian_previous_q_trunc, previous_q_trunc);
  gsl_vector_memcpy(previous_Gamma_trunc, &Gamma_trunc.vector);
  if (write_trajectory)
  {
    write_trajectory_output_files(condition, nb_iterations, t, dt);
  }
  while (nb_iterations < max_iter)
  {
    nb_iterations++;
    if (_stable_count >= stable_count)
    {
      break;
    }
    previous_mu = _mu;
    block_reactions();
    gsl_vector_view Gamma_trunc = gsl_vector_subvector(_Gamma, 1, _nj-1);
    /*
    if (hessian)
    {
      double alpha = 1.0;
      for (int j = 0; j < _nj-1; j++)
      {
        double gamma_j          = gsl_vector_get(&Gamma_trunc.vector, j);
        double previous_gamma_j = gsl_vector_get(previous_Gamma_trunc, j);
        double f_j              = gsl_vector_get(_q_trunc, j);
        double previous_f_j     = gsl_vector_get(hessian_previous_q_trunc, j);
        double gamma_diff       = fabs(gamma_j-previous_gamma_j);
        double f_diff           = fabs(f_j-previous_f_j);
        double h_j              = gamma_diff/f_diff;
        if (h_j < H_MIN)
        {
          h_j = H_MIN;
        }
        double h_j_inv = 1.0/h_j;
        if (h_j_inv > 1e+2)
        {
          h_j_inv = 1e+2;
        }
        double rescaled_gamma_j = gamma_j*h_j_inv;
        if (f_diff > H_MIN)
        {
          gsl_vector_set(scaled_Gammadt_trunc, j, alpha*rescaled_gamma_j+(1-alpha)*gamma_j);
        }
        else
        {
          gsl_vector_set(scaled_Gammadt_trunc, j, gamma_j);
        }
      }
    }
    else
    {
      gsl_vector_memcpy(scaled_Gammadt_trunc, &Gamma_trunc.vector);
    }
     */
    gsl_vector_memcpy(scaled_Gammadt_trunc, &Gamma_trunc.vector);
    gsl_vector_scale(scaled_Gammadt_trunc, dt);
    gsl_vector_add(_q_trunc, scaled_Gammadt_trunc);
    calculate_q_from_q_trunc();
    calculate();
    if (_consistent && _mu >= previous_mu)
    {
      gsl_vector_memcpy(hessian_previous_q_trunc, previous_q_trunc);
      gsl_vector_memcpy(previous_q_trunc, _q_trunc);
      gsl_vector_memcpy(previous_Gamma_trunc, &Gamma_trunc.vector);
      dt_counter++;
      t            = t+dt;
      _mu_diff     = fabs(_mu-previous_mu);
      _mu_rel_diff = fabs(_mu-previous_mu)/previous_mu;
      if (write_trajectory && nb_iterations%EXPORT_DATA_COUNT == 0)
      {
        save_q(nb_iterations, t, dt, output_path, condition);
        write_trajectory_output_files(condition, nb_iterations, t, dt);
        if (extra_verbose)
        {
          std::cout << " > Growth rate = " << _mu << " (iter=" << nb_iterations << ", mu_diff=" << _mu_diff << ", rel_diff=" << _mu_rel_diff << ", stable=" << _stable_count << ", dt=" << dt << ")" << std::endl;
        }
      }
      if (_mu_rel_diff < _mu_tol)
      {
        _stable_count++;
      }
      else
      {
        _stable_count--;
        if (_stable_count < 0)
        {
          _stable_count = 0;
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
      gsl_vector_memcpy(_q_trunc, previous_q_trunc);
      calculate_q_from_q_trunc();
      calculate();
      assert(_consistent);
      dt         /= DECREASING_DT_FACTOR;
      dt_counter  = 0;
      if (dt < 1e-100)
      {
        throw std::runtime_error("> Error: The timestep is too small (1e-100)");
      }
    }
  }
  save_q(nb_iterations, t, dt, output_path, condition);
  gsl_vector_free(previous_q_trunc);
  gsl_vector_free(previous_Gamma_trunc);
  gsl_vector_free(scaled_Gammadt_trunc);
  previous_q_trunc     = NULL;
  previous_Gamma_trunc = NULL;
  scaled_Gammadt_trunc = NULL;
  if (write_trajectory)
  {
    write_trajectory_output_files(condition, nb_iterations, t, dt);
    close_trajectory_ouput_files();
  }
  if (_stable_count >= stable_count)
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
 * \param    bool reload
 * \return   \e void
 */
void Model::open_trajectory_output_files( std::string output_path, std::string condition, bool reload )
{
  /*~~~~~~~~~~~~~~~~~~*/
  /* 1) Open files    */
  /*~~~~~~~~~~~~~~~~~~*/
  std::stringstream state_trajectory_filename;
  std::stringstream q_trajectory_filename;
  std::stringstream c_trajectory_filename;
  std::stringstream v_trajectory_filename;
  std::stringstream p_trajectory_filename;
  std::stringstream b_trajectory_filename;
  state_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_state_trajectory.csv";
  q_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_q_trajectory.csv";
  c_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_c_trajectory.csv";
  v_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_v_trajectory.csv";
  p_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_p_trajectory.csv";
  b_trajectory_filename << output_path << "/" << _model_name << "_" <<  condition << "_b_trajectory.csv";
  if (reload)
  {
    _state_trajectory_file.open(state_trajectory_filename.str(), std::ios::out | std::ios::app);
    _q_trajectory_file.open(q_trajectory_filename.str(), std::ios::out | std::ios::app);
    _c_trajectory_file.open(c_trajectory_filename.str(), std::ios::out | std::ios::app);
    _v_trajectory_file.open(v_trajectory_filename.str(), std::ios::out | std::ios::app);
    _p_trajectory_file.open(p_trajectory_filename.str(), std::ios::out | std::ios::app);
    _b_trajectory_file.open(b_trajectory_filename.str(), std::ios::out | std::ios::app);
  }
  else
  {
    _state_trajectory_file.open(state_trajectory_filename.str(), std::ios::out | std::ios::trunc);
    _q_trajectory_file.open(q_trajectory_filename.str(), std::ios::out | std::ios::trunc);
    _c_trajectory_file.open(c_trajectory_filename.str(), std::ios::out | std::ios::trunc);
    _v_trajectory_file.open(v_trajectory_filename.str(), std::ios::out | std::ios::trunc);
    _p_trajectory_file.open(p_trajectory_filename.str(), std::ios::out | std::ios::trunc);
    _b_trajectory_file.open(b_trajectory_filename.str(), std::ios::out | std::ios::trunc);
    /*~~~~~~~~~~~~~~~~~~*/
    /* 2) Write headers */
    /*~~~~~~~~~~~~~~~~~~*/
    _state_trajectory_file << "condition;iter;t;dt;mu;density;consistent;mu_diff;mu_rel_diff;stable_count\n";
    _q_trajectory_file << "condition;iter;t;dt";
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
      _q_trajectory_file << ";" << _reaction_ids[j];
      _v_trajectory_file << ";" << _reaction_ids[j];
      _p_trajectory_file << ";" << _reaction_ids[j];
    }
    _q_trajectory_file << "\n";
    _c_trajectory_file << "\n";
    _v_trajectory_file << "\n";
    _p_trajectory_file << "\n";
    _b_trajectory_file << "\n";
  }
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
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Update state file               */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  _state_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt << ";" << _mu << ";" << _density << ";" << _consistent << ";" << _mu_diff << ";" << _mu_rel_diff << ";" << _stable_count << "\n";
  _state_trajectory_file.flush();
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Update metabolites related file */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
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
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Update reactions related file   */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  _q_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt;
  _v_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt;
  _p_trajectory_file << condition << ";" << iter << ";" << t << ";" << dt;
  for (int j = 0; j < _nj; j++)
  {
    _q_trajectory_file << ";" << gsl_vector_get(_q, j);
    _v_trajectory_file << ";" << gsl_vector_get(_v, j);
    _p_trajectory_file << ";" << gsl_vector_get(_p, j);
  }
  _q_trajectory_file << "\n";
  _v_trajectory_file << "\n";
  _p_trajectory_file << "\n";
  _q_trajectory_file.flush();
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
  _q_trajectory_file.close();
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
  std::stringstream q_optimum_filename;
  std::stringstream c_optimum_filename;
  std::stringstream v_optimum_filename;
  std::stringstream p_optimum_filename;
  std::stringstream b_optimum_filename;
  state_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_state_optimum.csv";
  q_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_q_optimum.csv";
  c_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_c_optimum.csv";
  v_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_v_optimum.csv";
  p_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_p_optimum.csv";
  b_optimum_filename << output_path << "/" << _model_name << "_" << condition << "_b_optimum.csv";
  /*~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Create files     */
  /*~~~~~~~~~~~~~~~~~~~~~*/
  /*** Open files ***/
  _state_optimum_file.open(state_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _q_optimum_file.open(q_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _c_optimum_file.open(c_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _v_optimum_file.open(v_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _p_optimum_file.open(p_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _b_optimum_file.open(b_optimum_filename.str(), std::ios::out | std::ios::trunc);
  /*** Write headers ***/
  _state_optimum_file << "condition;mu;density;consistent;converged;run_time\n";
  _q_optimum_file << "condition";
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
    _q_optimum_file << ";" << _reaction_ids[j];
    _v_optimum_file << ";" << _reaction_ids[j];
    _p_optimum_file << ";" << _reaction_ids[j];
  }
  _q_optimum_file << "\n";
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
  _state_optimum_file << condition << ";" << _mu << ";" << _density << ";" << _consistent << ";" << converged << ";" << runtime << "\n";
  _q_optimum_file << condition;
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
    _q_optimum_file << ";" << gsl_vector_get(_q, j);
    _v_optimum_file << ";" << gsl_vector_get(_v, j);
    _p_optimum_file << ";" << gsl_vector_get(_p, j);
  }
  _q_optimum_file << "\n";
  _c_optimum_file << "\n";
  _v_optimum_file << "\n";
  _p_optimum_file << "\n";
  _b_optimum_file << "\n";
  _state_optimum_file.flush();
  _q_optimum_file.flush();
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
  _q_optimum_file.close();
  _c_optimum_file.close();
  _v_optimum_file.close();
  _p_optimum_file.close();
  _b_optimum_file.close();
}

/**
 * \brief    Save q vector into binary f
 * \details  --
 * \param    int nb_iterations
 * \param    double t
 * \param    double dt
 * \param    std::string output_path
 * \param    std::string condition
 * \return   \e void
 */
void Model::save_q( int nb_iterations, double t, double dt, std::string output_path, std::string condition )
{
  std::stringstream filename;
  filename << output_path << "/" << _model_name << "_" << condition << "_q.bin";
  FILE *f    = fopen(filename.str().c_str(), "wb");
  fwrite(&nb_iterations, sizeof(int), 1, f);
  fwrite(&t, sizeof(double), 1, f);
  fwrite(&dt, sizeof(double), 1, f);
  gsl_vector_fwrite(f, _q);
  fclose(f);
}

/**
 * \brief    Re-load q0
 * \details  Re-load the vector from the last trajectory point
 * \param    std::string output_path
 * \param    std::string condition
 * \param    bool restart
 * \return   \e void
 */
void Model::reload_q0( int &nb_iterations, double &t, double &dt, std::string output_path, std::string condition, bool restart )
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Check that q0 is empty           */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  gsl_vector_free(_q0);
  _q0 = gsl_vector_alloc(_nj);
  gsl_vector_set_zero(_q0);
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Load the binary file             */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  std::stringstream filename;
  filename << output_path << "/" << _model_name << "_" << condition << "_q.bin";
  FILE *f = fopen(filename.str().c_str(), "rb");
  fread(&nb_iterations, sizeof(int), 1, f);
  fread(&t, sizeof(double), 1, f);
  fread(&dt, sizeof(double), 1, f);
  gsl_vector_fread(f, _q0);
  fclose(f);
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Reset time variables if required */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  if (restart)
  {
    nb_iterations = 0;
    t             = 0.0;
    dt            = 0.01;
  }
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
  std::cout << "mu\tdensity\tconsistent\tconverged\trun_time" << std::endl;
  std::cout << _mu << "\t" << _density << "\t" << _consistent << "\t" << converged << "\t" << runtime << "\n";
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 3) f vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "variable";
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << _reaction_ids[j];
  }
  std::cout << std::endl << "q";
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << gsl_vector_get(_q, j);
  }
  std::cout << std::endl;
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 4) v vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "variable";
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << _reaction_ids[j];
  }
  std::cout << std::endl << "v";
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << gsl_vector_get(_v, j);
  }
  std::cout << std::endl;
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 5) p vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "variable";
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << _reaction_ids[j];
  }
  std::cout << std::endl << "p";
  for (int j = 0; j < _nj; j++)
  {
    std::cout << "\t" << gsl_vector_get(_p, j);
  }
  std::cout << std::endl;
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 6) b vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "variable";
  for (int j = 0; j < _nc; j++)
  {
    std::cout << "\t" << _c_ids[j];
  }
  std::cout << std::endl << "b";
  for (int j = 0; j < _nc; j++)
  {
    std::cout << "\t" << gsl_vector_get(_b, j);
  }
  std::cout << std::endl;
  /*~~~~~~~~~~~~~~~~~~~~*/
  /* 7) c vector output */
  /*~~~~~~~~~~~~~~~~~~~~*/
  std::cout << "variable";
  for (int j = 0; j < _nc; j++)
  {
    std::cout << "\t" << _c_ids[j];
  }
  std::cout << std::endl << "c";
  for (int j = 0; j < _nc; j++)
  {
    std::cout << "\t" << gsl_vector_get(_c, j);
  }
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
    id.erase(std::remove(id.begin(), id.end(), '\r'), id.end());
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
    id.erase(std::remove(id.begin(), id.end(), '\r'), id.end());
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
 * \brief    Load the K matrix
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_K( void )
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Load the complete KM matrix             */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  assert(_K==NULL);
  _K = gsl_matrix_alloc(_ni, _nj);
  gsl_matrix_set_zero(_K);
  int row = 0;
  int col = 0;
  assert(is_file_exist(_model_path+"/"+_model_name+"/K.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/K.csv", std::ios::in);
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
      gsl_matrix_set(_K, row, col, value);
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
        gsl_matrix_set(_KM_f, i, j, gsl_matrix_get(_K, i, j));
      }
      else if (gsl_matrix_get(_Mx, i, j) > 0.0)
      {
        gsl_matrix_set(_KM_b, i, j, gsl_matrix_get(_K, i, j));
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
    id.erase(std::remove(id.begin(), id.end(), '\r'), id.end());
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
    id.erase(std::remove(id.begin(), id.end(), '\r'), id.end());
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
  if (is_file_exist(_model_path+"/"+_model_name+"/constant_reactions.csv"))
  {
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
      reaction_id.erase(std::remove(reaction_id.begin(), reaction_id.end(), '\r'), reaction_id.end());
      assert(_constant_reactions.find(reaction_id) == _constant_reactions.end());
      _constant_reactions[reaction_id] = stod(str_value);
    }
    file.close();
  }
}

/**
 * \brief    Load q0
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_q0( void )
{
  assert(_q0==NULL);
  _q0 = gsl_vector_alloc(_nj);
  gsl_vector_set_zero(_q0);
  assert(is_file_exist(_model_path+"/"+_model_name+"/q.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/q.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string str_value;
  getline(file, line);
  getline(file, line);
  std::stringstream flux(line.c_str());
  int index = -1;
  while(getline(flux, str_value, ';'))
  {
    if (index==-1)
    {
      assert(str_value=="q0");
    }
    else
    {
      double value = stod(str_value);
      gsl_vector_set(_q0, index, value);
    }
    index++;
  }
  file.close();
  assert(index==_nj);
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
  for(int j = 0; j < _nj; j++)
  {
    gsl_matrix_get_col(KI_vec, _KI, j);
    gsl_matrix_get_col(KA_vec, _KA, j);
    bool kcat_b_zero = fabs(gsl_vector_get(_kcat_b, j)) < _tol;
    bool KI_sum_zero = fabs(gsl_blas_dasum(KI_vec)) < _tol;
    bool KA_sum_zero = fabs(gsl_blas_dasum(KA_vec)) < _tol;
    if (kcat_b_zero && KI_sum_zero && KA_sum_zero)
    {
      _type[j] = IMM;
    }
    else if (kcat_b_zero && !KI_sum_zero && KA_sum_zero)
    {
      _type[j] = IMMI;
    }
    else if (kcat_b_zero && KI_sum_zero && !KA_sum_zero)
    {
      _type[j] = IMMA;
    }
    else if (kcat_b_zero && !KI_sum_zero && !KA_sum_zero)
    {
      _type[j] = IMMIA;
    }
    else if (!kcat_b_zero)
    {
      assert(KI_sum_zero);
      assert(KA_sum_zero);
      _type[j] = RMM;
    }
    else
    {
      throw std::runtime_error("> Error: The kinetic scheme of reaction "+_reaction_ids[j]+" is incorrect");
    }
  }
  gsl_vector_free(KI_vec);
  gsl_vector_free(KA_vec);
  KI_vec = NULL;
  KA_vec = NULL;
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
  assert(_q==NULL);
  assert(_q_trunc==NULL);
  assert(_c==NULL);
  assert(_xc==NULL);
  assert(_tau_j==NULL);
  assert(_v==NULL);
  assert(_p==NULL);
  assert(_b==NULL);
  assert(_ditau_j==NULL);
  assert(_dmu_dq==NULL);
  assert(_Gamma==NULL);
  assert(_dmu_dq_term2==NULL);
  assert(_dmu_dq_term3==NULL);
  assert(_dmu_dq_term4==NULL);
  assert(_dmu_dq_term5==NULL);
  /*** Allocate memory ***/
  _x            = gsl_vector_alloc(_nx);
  _q            = gsl_vector_alloc(_nj);
  _q_trunc      = gsl_vector_alloc(_nj-1);
  _c            = gsl_vector_alloc(_nc);
  _xc           = gsl_vector_alloc(_ni);
  _tau_j        = gsl_vector_alloc(_nj);
  _v            = gsl_vector_alloc(_nj);
  _p            = gsl_vector_alloc(_nj);
  _b            = gsl_vector_alloc(_nc);
  _ditau_j      = gsl_matrix_alloc(_nj, _nc);
  _dmu_dq       = gsl_vector_alloc(_nj);
  _Gamma        = gsl_vector_alloc(_nj);
  _dmu_dq_term2 = gsl_vector_alloc(_nj);
  _dmu_dq_term3 = gsl_matrix_alloc(_nj, _nj);
  _dmu_dq_term4 = gsl_vector_alloc(_nj);
  _dmu_dq_term5 = gsl_vector_alloc(_nj);
  /*** Initialize all variables to zero ***/
  gsl_vector_set_zero(_x);
  gsl_vector_set_zero(_q);
  gsl_vector_set_zero(_q_trunc);
  gsl_vector_set_zero(_c);
  gsl_vector_set_zero(_xc);
  gsl_vector_set_zero(_tau_j);
  gsl_vector_set_zero(_v);
  gsl_vector_set_zero(_p);
  gsl_vector_set_zero(_b);
  gsl_matrix_set_zero(_ditau_j);
  gsl_vector_set_zero(_dmu_dq);
  gsl_vector_set_zero(_Gamma);
  gsl_vector_set_zero(_dmu_dq_term2);
  gsl_matrix_set_zero(_dmu_dq_term3);
  gsl_vector_set_zero(_dmu_dq_term4);
  gsl_vector_set_zero(_dmu_dq_term5);
  /*** Initialize vector views ***/
  _x_view = gsl_vector_subvector(_xc, 0, _nx);
  _c_view = gsl_vector_subvector(_xc, _nx, _nc);
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
  gsl_blas_dgemv(CblasNoTrans, _rho, _M, _q, 0.0, _c);
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
  // iMM <- function(j,c,xc) as.numeric(prod(1 + KS[,j]/xc)/kcatf[j])
  double prod_KM_f = 1.0;
  double kcatf     = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    prod_KM_f *= 1.0 + gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  gsl_vector_set(_tau_j, j, prod_KM_f/kcatf);
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + inhibition (only one inhibitor per reaction)
 * \details  Formula: tau_j = prod(1+xc*1/KI[,j])*prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMMi( int j )
{
  // iMMi <- function(j,c,xc) as.numeric( prod(1 + xc*rKI[,j])*prod(1 + KS[,j]/xc) /kcatf[j] )
  double prod_KI   = 1.0;
  double prod_KM_f = 1.0;
  double kcatf     = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    double rKI = (gsl_matrix_get(_KI, i, j) > _tol ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    prod_KI   *= 1.0 + gsl_vector_get(_xc, i)*rKI;
    prod_KM_f *= 1.0 + gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  gsl_vector_set(_tau_j, j, prod_KI*prod_KM_f/kcatf);
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + activation (only one activator per reaction)
 * \details  Formula: tau_j = prod(1+KA[,j]/xc)*prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMMa( int j )
{
  // iMMa <- function(j,c,xc) as.numeric(prod(1 + KA[,j]/xc)*prod(1 + KS[,j]/xc)/kcatf[j])
  double prod_KA   = 1.0;
  double prod_KM_f = 1.0;
  double kcatf     = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    prod_KA   *= 1.0 + gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
    prod_KM_f *= 1.0 + gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  
  gsl_vector_set(_tau_j, j, prod_KA*prod_KM_f/kcatf);
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + inhibition + activation
 * \details  Formula: tau_j = prod(1+xc*1/KI[,j])*prod(1+KA[,j]/xc)*prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMMia( int j )
{
  // iMMia <- function(j,c,xc)  as.numeric(prod(1 + xc*rKI[,j])*prod(1 + KA[,j]/xc)*prod(1 + KS[,j]/xc)/kcat[j])
  double prod_KI   = 1.0;
  double prod_KA   = 1.0;
  double prod_KM_f = 1.0;
  double kcatf     = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    double rKI  = (gsl_matrix_get(_KI, i, j) > _tol ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    prod_KI    *= 1.0 + gsl_vector_get(_xc, i)*rKI;
    prod_KA    *= 1.0 + gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
    prod_KM_f  *= 1.0 + gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  gsl_vector_set(_tau_j, j, prod_KI*prod_KA*prod_KM_f/kcatf);
}

/**
 * \brief    Reversible Michaelis-Menten kinetics
 * \details  Formula: tau_j = 1/[ kcat_f[j]/prod(1+Km_f[,j]/xc) - kcat_b[j]/prod(1+Km_b[,j]/xc) ]
 * \param    int j
 * \return   \e void
 */
void Model::rMM( int j )
{
  // rMM <- 1 / ( kcatf[j]/prod(1 + KS[,j]/xc) - kcatb[j]/prod(1 + KP[,j]/xc)  )
  double prod_KM_f = 1.0;
  double prod_KM_b = 1.0;
  double kcatf     = gsl_vector_get(_kcat_f, j);
  double kcatb     = gsl_vector_get(_kcat_b, j);
  for (int i = 0; i < _ni; i++)
  {
    prod_KM_f *= 1.0 + gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    prod_KM_b *= 1.0 + gsl_matrix_get(_KM_b, i, j)/gsl_vector_get(_xc, i);
  }
  gsl_vector_set(_tau_j, j, 1.0/(kcatf/prod_KM_f-kcatb/prod_KM_b));
}

/**
 * \brief    Global Michaelis-Menten kinetics
 * \details  Formula: tauj[j] = prod(1 + KA[,j]/ca)*prod(1 + ca/KI[,j])*( prod(1 + subr) + prod(1 + prodr) - 1 )/( kcatf[j]*prod(subr) - kcatb[j]*prod(prodr)  )
 * \param    int j
 * \return   \e void
 */
/*
void Model::gMM( int j )
{
  double kcatf       = gsl_vector_get(_kcat_f, j);
  double kcatb       = gsl_vector_get(_kcat_b, j);
  double prod_KM_f   = 1.0;
  double prod_KM_b   = 1.0;
  double prod_KM_f_1 = 1.0;
  double prod_KM_b_1 = 1.0;
  double prod_KA_1   = 1.0;
  double prod_KI_1   = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    double x_c  = gsl_vector_get(_xc, i);
    double KM_f = gsl_matrix_get(_KM_f, i, j);
    double KM_b = gsl_matrix_get(_KM_b, i, j);
    double KA   = gsl_matrix_get(_KA, i, j);
    double rKI  = (gsl_matrix_get(_KI, i, j) > _tol ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    if (KM_f > 0.0)
    {
      prod_KM_f   *= x_c/KM_f;
      prod_KM_f_1 *= 1.0 + x_c/KM_f;
    }
    if (KM_b > 0.0)
    {
      prod_KM_b   *= x_c/KM_b;
      prod_KM_b_1 *= 1.0 + x_c/KM_b;
    }
    prod_KA_1 *= 1.0 + KA/x_c;
    prod_KI_1 *= 1.0 + x_c*rKI;
  }
  double tau_j = prod_KA_1*prod_KI_1*(prod_KM_f_1+prod_KM_b_1-1.0)/(kcatf*prod_KM_f-kcatb*prod_KM_b);
  gsl_vector_set(_tau_j, j, tau_j);
}
*/

/**
 * \brief    Compute tau_j
 * \details  --
 * \param    int j
 * \return   \e void
 */
void Model::compute_tau( int j )
{
  //gMM(j);
  //return;
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
  double kcatf = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double term1 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term2 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term2 *= 1.0 + gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
    }
    gsl_matrix_set(_ditau_j, j, i, -term1*term2/kcatf);
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
  double prod_KI   = 1.0;
  double prod_KM_f = 1.0;
  double kcatf     = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    double rKI  = (gsl_matrix_get(_KI, i, j) > _tol ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    prod_KI    *= 1.0 + gsl_vector_get(_xc, i)*rKI;
    prod_KM_f  *= 1.0 + gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    
  }
  // ditauj[i2] <- ( rKI[y,j] * prod_KM_f - prod_KI * (KS[y,j]/(c[i2]^2)) * prod(1 + KS[-y,j]/xc[-y]) )/kcatf[j]
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double rKI   = (gsl_matrix_get(_KI, y, j) > _tol ? 1.0/gsl_matrix_get(_KI, y, j) : 0.0);
    double term1 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term2 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term2 *= 1.0 + gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
    }
    gsl_matrix_set(_ditau_j, j, i, (rKI*prod_KM_f-prod_KI*term1*term2)/kcatf);
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
  double prod_KA   = 1.0;
  double prod_KM_f = 1.0;
  double kcatf     = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    prod_KA   *= 1.0 + gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
    prod_KM_f *= 1.0 + gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  // ditauj[i2] <- -as.numeric( term1*prod_KM_f + term2*prod_KA*term3 )/kcatf[j]
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
        term3 *= 1.0 + gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
    }
    gsl_matrix_set(_ditau_j, j, i, -(term1*prod_KM_f+term2*prod_KA*term3)/kcatf);
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
  double prod_KI   = 1.0;
  double prod_KA   = 1.0;
  double prod_KM_f = 1.0;
  double kcatf     = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    double rKI  = (gsl_matrix_get(_KI, i, j) > _tol ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    prod_KI    *= 1.0 + gsl_vector_get(_xc, i)*rKI;
    prod_KA    *= 1.0 + gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
    prod_KM_f  *= 1.0 + gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double rKI   = (gsl_matrix_get(_KI, y, j) > _tol ? 1.0/gsl_matrix_get(_KI, y, j) : 0.0);
    double term2 = -gsl_matrix_get(_KA, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term3 = -gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term4 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term4 *= 1.0 + gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
    }
    gsl_matrix_set(_ditau_j, j, i, (rKI*prod_KA*prod_KM_f+prod_KI*term2*prod_KM_f+prod_KI*prod_KA*term3*term4)/kcatf);
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
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Calculate constant terms for reaction j */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  double prod_KM_f = 1.0;
  double prod_KM_b = 1.0;
  double kcatf     = gsl_vector_get(_kcat_f, j);
  double kcatb     = gsl_vector_get(_kcat_b, j);
  for (int i = 0; i < _ni; i++)
  {
    prod_KM_f *= 1.0 + gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    prod_KM_b *= 1.0 + gsl_matrix_get(_KM_b, i, j)/gsl_vector_get(_xc, i);
  }
  //double tau_j_2 = 1.0 / gsl_pow_int(kcatf/prod_KM_f-kcatb/prod_KM_b, 2);
  double tau_j = 1.0/(kcatf/prod_KM_f-kcatb/prod_KM_b);
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Calculate terms depending on substrate  */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  for (int i = 0; i < _nc; i++)
  {
    int    y      = i+_nx;
    double term1 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i)+gsl_matrix_get(_KM_f, y, j), 2);
    double term2 = gsl_matrix_get(_KM_b, y, j)/gsl_pow_int(gsl_vector_get(_c, i)+gsl_matrix_get(_KM_b, y, j), 2);
    double prodf  = 1.0;
    double prodb  = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        prodf *= 1.0 + gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
        prodb *= 1.0 + gsl_matrix_get(_KM_b, index, j)/gsl_vector_get(_xc, index);
      }
    }
    double ditauj = (kcatf/prodf)*term1 - (kcatb/prodb)*term2;
    gsl_matrix_set(_ditau_j, j, i, -ditauj*tau_j*tau_j);
  }
}

/**
 * \brief    Global derivative function with respect to metabolite concentrations
 * \details  Formula: --
 * \param    int j
 * \return   \e void
 */
/*
void Model::dgMM( int j )
{
  double kcatf       = gsl_vector_get(_kcat_f, j);
  double kcatb       = gsl_vector_get(_kcat_b, j);
  double prod_KM_f   = 1.0;
  double prod_KM_b   = 1.0;
  double prod_KM_f_1 = 1.0;
  double prod_KM_b_1 = 1.0;
  double prod_KA_1   = 1.0;
  double prod_KI_1   = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    double x_c  = gsl_vector_get(_xc, i);
    double KM_f = gsl_matrix_get(_KM_f, i, j);
    double KM_b = gsl_matrix_get(_KM_b, i, j);
    double KA   = gsl_matrix_get(_KA, i, j);
    double rKI  = (gsl_matrix_get(_KI, i, j) > 0.0 ? 1.0/gsl_matrix_get(_KI, i, j) : 0.0);
    if (KM_f > 0.0)
    {
      prod_KM_f   *= x_c/KM_f;
      prod_KM_f_1 *= 1.0 + x_c/KM_f;
    }
    if (KM_b > 0.0)
    {
      prod_KM_b   *= x_c/KM_b;
      prod_KM_b_1 *= 1.0 + x_c/KM_b;
    }
    prod_KA_1 *= 1.0 + KA/x_c;
    prod_KI_1 *= 1.0 + x_c*rKI;
  }
  double tau_j_A = prod_KA_1*prod_KI_1;
  double tau_j_B = prod_KM_f_1+prod_KM_b_1-1.0;
  double tau_j_C = kcatf*prod_KM_f-kcatb*prod_KM_b;
  for (int i = 0; i < _nc; i++)
  {
    int y = i+_nx;
    //### 2.1) Calculate products depending on y ###
    double prod_KM_f_y   = 1.0;
    double prod_KM_b_y   = 1.0;
    double prod_KM_f_1_y = 1.0;
    double prod_KM_b_1_y = 1.0;
    double prod_KA_1_y   = 1.0;
    double prod_KI_1_y   = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        double rKI  = (gsl_matrix_get(_KI, index, j) > 0.0 ? 1.0/gsl_matrix_get(_KI, index, j) : 0.0);
        double x_c  = gsl_vector_get(_xc, index);
        double KM_f = gsl_matrix_get(_KM_f, index, j);
        double KM_b = gsl_matrix_get(_KM_b, index, j);
        if (KM_f > 0.0)
        {
          prod_KM_f_y   *= x_c/KM_f;
          prod_KM_f_1_y *= 1.0 + x_c/KM_f;
        }
        if (KM_b > 0.0)
        {
          prod_KM_b_y   *= x_c/KM_b;
          prod_KM_b_1_y *= 1.0 + x_c/KM_b;
        }
        prod_KA_1_y *= 1.0 + gsl_matrix_get(_KA, index, j)/x_c;
        prod_KI_1_y *= 1.0 + x_c*rKI;
      }
    }
    //### 2.2) Calculate dtau_j terms ###
    double KM_f  = gsl_matrix_get(_KM_f, y, j);
    double KM_b  = gsl_matrix_get(_KM_b, y, j);
    double rKI   = (gsl_matrix_get(_KI, y, j) > 0.0 ? 1.0/gsl_matrix_get(_KI, y, j) : 0.0);
    double term1 = gsl_matrix_get(_KA, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term2 = 0.0;
    double term3 = 0.0;
    double term4 = 0.0;
    double term5 = 0.0;
    if (KM_f > 0.0)
    {
      term2 = prod_KM_f_1_y/KM_f;
      term4 = kcatf/KM_f*prod_KM_f_y;
    }
    if (KM_b > 0.0)
    {
      term3 = prod_KM_b_1_y/KM_b;
      term5 = kcatb/KM_b*prod_KM_b_y;
    }
    double dtau_j_A = -term1*prod_KA_1_y*prod_KI_1 + prod_KA_1*rKI*prod_KI_1_y;
    double dtau_j_B = term2+term3;
    double dtau_j_C = term4-term5;
    //### 2.3) Calculate dtau_j ###
    double ditau_j = ( dtau_j_A*tau_j_B + tau_j_A*dtau_j_B - tau_j_A*tau_j_B*dtau_j_C/tau_j_C  )/tau_j_C;
    gsl_matrix_set(_ditau_j, j, i, ditau_j);
  }
}
*/

/**
 * \brief    Compute ditau_j
 * \details  --
 * \param    int j
 * \return   \e void
 */
void Model::compute_dtau( int j )
{
  //dgMM(j);
  //return;
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
  gsl_blas_ddot(_tau_j, _q, &term1);
  _mu = gsl_matrix_get(_M, _a, _r)*gsl_vector_get(_q, _r)/term1;
}

/**
 * \brief    Compute the mass flux vector v
 * \details  v = mu*rho*q
 * \param    void
 * \return   \e void
 */
void Model::compute_v( void )
{
  gsl_vector_memcpy(_v, _q);
  gsl_vector_scale(_v, _mu*_rho);
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
 * \details  b = M*q
 * \param    void
 * \return   \e void
 */
void Model::compute_b( void )
{
  gsl_blas_dgemv(CblasNoTrans, 1.0, _M, _q, 0.0, _b);
}

/**
 * \brief    Compute cell density (should be always 1)
 * \details  density = sM*q
 * \param    void
 * \return   \e void
 */
void Model::compute_density( void )
{
  gsl_blas_ddot(_sM, _q, &_density);
}

/**
 * \brief    Compute local mu gradient with respect to q
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::compute_dmu_dq( void )
{
  // ((mu(f)^2)/b(f)[p]) * (M[p,]/mu(f)   - t(f) %*% ( rho*dtau(ci(f))%*%M ) - tau(ci(f)) )
  /*--------*/
  _dmu_dq_term1 = gsl_pow_int(_mu, 2)/gsl_vector_get(_b, _a);
  /*--------*/
  gsl_matrix_get_row(_dmu_dq_term2, _M, _a);
  gsl_vector_scale(_dmu_dq_term2, 1.0/_mu);
  /*--------*/
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, _rho, _ditau_j, _M, 0.0, _dmu_dq_term3);
  gsl_blas_dgemv(CblasTrans, 1.0, _dmu_dq_term3, _q, 0.0, _dmu_dq_term4);
  /*--------*/
  gsl_vector_memcpy(_dmu_dq_term5, _tau_j);
  /*--------*/
  gsl_vector_memcpy(_dmu_dq, _dmu_dq_term2);
  gsl_vector_sub(_dmu_dq, _dmu_dq_term4);
  gsl_vector_sub(_dmu_dq, _dmu_dq_term5);
  gsl_vector_scale(_dmu_dq, _dmu_dq_term1);
}

/**
 * \brief    Compute local growth control coefficients with respect to q
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::compute_Gamma( void )
{
  gsl_vector_memcpy(_Gamma, _sM);
  gsl_vector_scale(_Gamma, gsl_vector_get(_dmu_dq, 0)/gsl_vector_get(_sM, 0));
  gsl_vector_scale(_Gamma, -1.0);
  gsl_vector_add(_Gamma, _dmu_dq);
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
  compute_dmu_dq();
  compute_Gamma();
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
  if (gsl_vector_min(_c) < 0.0)
  {
    test2 = false;
  }
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Test negative proteins constraint       */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  bool test3 = true;
  if (gsl_vector_min(_p) < 0.0)
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
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Manage reactions converging to zero */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  for (int j = 0; j < _nj-1; j++)
  {
    /*** 1.1) Reaction is irreversible ***/
    if (gsl_vector_get(_q_trunc, j) <= _tol) // (_type[j+1] != RMM && gsl_vector_get(_q_trunc, j) <= _tol)
    {
      gsl_vector_set(_q_trunc, j, _tol);
      if (gsl_vector_get(_Gamma, j+1) < 0.0)
      {
        gsl_vector_set(_Gamma, j+1, 0.0);
      }
    }
    /*** 1.2) Reaction is reversible ***/
    /*
    else if (_type[j+1] == RMM && fabs(gsl_vector_get(_q_trunc, j)) < _tol)
    {
      double gamma_j = gsl_vector_get(_Gamma, j+1);
      double q_j     = gsl_vector_get(_q_trunc, j);
      std::cout << q_j << " " << gamma_j << "\n";
      if (q_j < 0.0 && q_j > -_tol && gamma_j > 0.0)
      {
        gsl_vector_set(_q_trunc, j, _tol);
      }
      if (q_j > 0.0 && q_j < _tol && gamma_j < 0.0)
      {
        gsl_vector_set(_q_trunc, j, -_tol);
      }
    }
     */
  }
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Manage constant reactions           */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  for (auto item : _constant_reactions)
  {
    int    j   = _reaction_indices[item.first];
    double val = item.second;
    assert(j > 0);
    gsl_vector_set(_q_trunc, j-1, val);
    gsl_vector_set(_Gamma, j, 0.0);
  }
}

