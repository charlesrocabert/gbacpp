/**
 * \file      Model.cpp
 * \author    Charles Rocabert
 * \date      22-07-2024
 * \copyright GBA_Evolution. Copyright © 2024 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Model class definition
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

#include "Model.hpp"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    std::string model_path
 * \param    std::string model_name
 * \param    msize model_size
 * \return   \e void
 */
Model::Model( std::string model_path, std::string model_name, msize model_size )
{
  /*----------------------------------------------- Model path and name */
  
  _model_path = model_path;
  _model_name = model_name;
  _model_size = model_size;
  
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
  _KM_f       = NULL;
  _KM_b       = NULL;
  _KI         = NULL;
  _KA         = NULL;
  _kcat_f     = NULL;
  _kcat_b     = NULL;
  _type       = NULL;
  _conditions = NULL;
  _directions = NULL;
  _boundaries = NULL;

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
  
  _f0         = NULL;
  _f          = NULL;
  _f_trunc    = NULL;
  _c          = NULL;
  _xc         = NULL;
  _tau_j      = NULL;
  _v          = NULL;
  _p          = NULL;
  _b          = NULL;
  _density    = 0.0;
  _mu         = 0.0;
  _consistent = false;
  
  /*----------------------------------------------- GBA second order variables */
  
  _ditau_j = NULL;
  _dmu_f   = NULL;
  _GCC_f   = NULL;
  
  /*----------------------------------------------- Variables for calculation optimization */
  
  _KM_f_product = 1.0;
  _KM_b_product = 1.0;
  _KI_product   = 1.0;
  _KA_product   = 1.0;
  _dmu_f_term1  = 0.0;
  _dmu_f_term2  = NULL;
  _dmu_f_term3  = NULL;
  _dmu_f_term4  = NULL;
  _dmu_f_term5  = NULL;
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
  gsl_matrix_free(_KM_f);
  gsl_matrix_free(_KM_b);
  gsl_matrix_free(_KI);
  gsl_matrix_free(_KA);
  gsl_vector_free(_kcat_f);
  gsl_vector_free(_kcat_b);
  delete[] _type;
  gsl_matrix_free(_conditions);
  delete[] _directions;
  delete[] _boundaries;
  _Mx         = NULL;
  _M          = NULL;
  _KM_f       = NULL;
  _KM_b       = NULL;
  _KI         = NULL;
  _KA         = NULL;
  _kcat_f     = NULL;
  _kcat_b     = NULL;
  _type       = NULL;
  _conditions = NULL;
  _directions = NULL;
  _boundaries = NULL;
  
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
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Load the model from CSV files
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_model( void )
{
  load_metabolite_identifiers();
  load_reaction_identifiers();
  load_vector_sizes();
  load_M();
  load_KM_forward();
  load_KM_backward();
  load_KI();
  load_KA();
  load_kcat();
  load_conditions();
  load_f0();
  load_directions_and_boundaries();
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
 * \brief    Calculate all variables for genome-scale models
 * \details  Specific methods are applied
 * \param    void
 * \return   \e void
 */
void Model::calculate_GM( void )
{
  compute_c();
  compute_xc();
  for (int j = 0; j < _nj; j++)
  {
    GMMM(j);
  }
  compute_mu();
  compute_v();
  compute_p();
  compute_b();
  compute_density();
  for (int j = 0; j < _nj; j++)
  {
    dGMMM(j);
  }
  compute_dmu_f();
  compute_GCC_f();
  check_model_consistency();
}

/**
 * \brief    Calculate first order variables for genome-scale models
 * \details  Specific methods are applied
 * \param    void
 * \return   \e void
 */
void Model::calculate_first_order_GM( void )
{
  compute_c();
  compute_xc();
  for (int j = 0; j < _nj; j++)
  {
    GMMM(j);
  }
  compute_mu();
  compute_v();
  compute_p();
  compute_b();
  compute_density();
  check_model_consistency();
}

/**
 * \brief    Calculate second order variables for genome-scale models
 * \details  Specific methods are applied
 * \param    void
 * \return   \e void
 */
void Model::calculate_second_order_GM( void )
{
  compute_c();
  compute_xc();
  for (int j = 0; j < _nj; j++)
  {
    dGMMM(j);
  }
  compute_dmu_f();
  compute_GCC_f();
}

/**
 * \brief    Compute the gradient ascent trajectory
 * \details  --
 * \param    std::string condition
 * \param    double initial_dt
 * \param    double max_t
 * \param    double save_trajectory
 * \return   \e bool
 */
bool Model::compute_gradient_ascent_trajectory( std::string condition, double initial_dt, double max_t, double save_trajectory )
{
  if (_model_size == SMALL)
  {
    return compute_gradient_ascent_trajectory_for_small_models(condition, initial_dt, max_t, save_trajectory);
  }
  else if (_model_size == GENOME_SCALE)
  {
    return compute_gradient_ascent_trajectory_for_genome_scale_models(condition, initial_dt, max_t, save_trajectory);
  }
  else
  {
    throw std::invalid_argument("> Model size parameter incorrect value");
  }
  return false;
}

/**
 * \brief    Compute local optimum for all conditions
 * \details  --
 * \param    std::string condition
 * \param    double initial_dt
 * \param    double max_t
 * \return   \e bool
 */
void Model::compute_local_optimum_for_all_conditions( double initial_dt, double max_t )
{
  open_optimum_output_files();
  for (int i = 0; i < _condition_ids.size(); i++)
  {
    std::string condition = _condition_ids[i];
    bool converged        = compute_gradient_ascent_trajectory(condition, initial_dt, max_t, false);
    write_optimum_output_files(condition, converged);
  }
  close_optimum_ouput_files();
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Test file existence
 * \details  --
 * \param    std::string filename
 * \return   \e bool
 */
bool Model::is_file_exist( std::string filename )
{
    std::ifstream infile(filename.c_str());
    return infile.good();
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
  line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
  int index = 0;
  while(getline(file, line))
  {
    line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
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
  line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
  std::stringstream flux(line.c_str());
  getline(flux, id, ';');
  int index = 0;
  while(getline(flux, id, ';'))
  {
    line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
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
    line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
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
 * \brief    Load the KM forward matrix
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_KM_forward( void )
{
  assert(_KM_f==NULL);
  _KM_f = gsl_matrix_alloc(_ni, _nj);
  gsl_matrix_set_zero(_KM_f);
  int row = 0;
  int col = 0;
  assert(is_file_exist(_model_path+"/"+_model_name+"/KM_forward.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/KM_forward.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string id;
  std::string str_value;
  /*** skip header line ***/
  getline(file, line);
  while(getline(file, line))
  {
    line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
    std::stringstream flux(line.c_str());
    /*** get metabolite id ***/
    getline(flux, id, ';');
    col = 0;
    /*** fill the matrices for the given line ***/
    while(getline(flux, str_value, ';'))
    {
      double value = stod(str_value);
      gsl_matrix_set(_KM_f, row, col, value);
      col++;
    }
    row++;
  }
  file.close();
}

/**
 * \brief    Load the KM backward matrix
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_KM_backward( void )
{
  assert(_KM_b==NULL);
  _KM_b = gsl_matrix_alloc(_ni, _nj);
  gsl_matrix_set_zero(_KM_b);
  if (is_file_exist(_model_path+"/"+_model_name+"/KM_backward.csv"))
  {
    int row = 0;
    int col = 0;
    assert(is_file_exist(_model_path+"/"+_model_name+"/KM_backward.csv"));
    std::ifstream file(_model_path+"/"+_model_name+"/KM_backward.csv", std::ios::in);
    assert(file);
    std::string line;
    std::string id;
    std::string str_value;
    /*** skip header line ***/
    getline(file, line);
    while(getline(file, line))
    {
      line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
      std::stringstream flux(line.c_str());
      /*** get metabolite id ***/
      getline(flux, id, ';');
      col = 0;
      /*** fill the matrices for the given line ***/
      while(getline(flux, str_value, ';'))
      {
        double value = stod(str_value);
        gsl_matrix_set(_KM_b, row, col, value);
        col++;
      }
      row++;
    }
    file.close();
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
      line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
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
      line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
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
  line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
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
  line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
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
  /*-------------------------------*/
  /* 1) Read condition identifiers */
  /*-------------------------------*/
  getline(file, line);
  line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
  std::stringstream flux1(line.c_str());
  getline(flux1, id, ';');
  int index = 0;
  while(getline(flux1, id, ';'))
  {
    _condition_ids.push_back(id);
    _condition_indices[id] = index;
    index++;
  }
  /*-------------------------------*/
  /* 2) Read condition parameters  */
  /*-------------------------------*/
  while(getline(file, line))
  {
    line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
    std::stringstream flux(line.c_str());
    getline(flux, id, ';');
    _condition_params.push_back(id);
  }
  file.close();
  /*-------------------------------*/
  /* 3) Load the conditions        */
  /*-------------------------------*/
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
    if (fabs(value) > FLUX_BOUNDARY)
    {
      throw std::invalid_argument("> f0 value is higher than the flux boundary");
    }
  }
  file.close();
}

/**
 * \brief    Load directions and boundaries
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::load_directions_and_boundaries( void )
{
  assert(_directions==NULL);
  assert(_boundaries==NULL);
  _directions = new std::string[_nj];
  _boundaries = new boundaries[_nj];
  assert(is_file_exist(_model_path+"/"+_model_name+"/direction.csv"));
  std::ifstream file(_model_path+"/"+_model_name+"/direction.csv", std::ios::in);
  assert(file);
  std::string line;
  std::string reaction_id;
  std::string direction;
  getline(file, line);
  while(getline(file, line))
  {
    std::stringstream flux(line.c_str());
    getline(flux, reaction_id, ';');
    getline(flux, direction, ';');
    assert(_reaction_indices.find(reaction_id) != _reaction_indices.end());
    _directions[_reaction_indices[reaction_id]] = direction;
    if (direction == "forward")
    {
      _boundaries[_reaction_indices[reaction_id]].LB = 0.0;
      _boundaries[_reaction_indices[reaction_id]].UB = FLUX_BOUNDARY;
    }
    else if (direction == "backward")
    {
      _boundaries[_reaction_indices[reaction_id]].LB = -FLUX_BOUNDARY;
      _boundaries[_reaction_indices[reaction_id]].UB = 0.0;
    }
    else if (direction == "reversible")
    {
      _boundaries[_reaction_indices[reaction_id]].LB = -FLUX_BOUNDARY;
      _boundaries[_reaction_indices[reaction_id]].UB = FLUX_BOUNDARY;
    }
    else
    {
      throw std::invalid_argument("> Incorrect direction value");
    }
  }
  file.close();
}

/**
 * \brief    Initialize static variables
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::initialize_static_variables( void )
{
  /*------------------------------*/
  /* 1) Initialize constants      */
  /*------------------------------*/
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
  /*------------------------------*/
  /* 2) Initialize reaction types */
  /*------------------------------*/
  _type              = new rtype[_nj];
  gsl_vector* KI_vec = gsl_vector_alloc(_ni);
  gsl_vector* KA_vec = gsl_vector_alloc(_ni);
  for(int j = 0; j < _nj; j++)
  {
    gsl_matrix_get_col(KI_vec, _KI, j);
    gsl_matrix_get_col(KA_vec, _KA, j);
    bool kcat_b_zero = fabs(gsl_vector_get(_kcat_b, j)) < MIN_PARAMETER;
    bool KI_sum_zero = fabs(gsl_blas_dasum(KI_vec)) < MIN_PARAMETER;
    bool KA_sum_zero = fabs(gsl_blas_dasum(KA_vec)) < MIN_PARAMETER;
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
      _type[j] = IMMIA;;
    }
    else if (!kcat_b_zero)
    {
      assert(KI_sum_zero);
      assert(KA_sum_zero);
      _type[j] = RMM;
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
  _dmu_f_term3 = gsl_matrix_alloc(_nc, _nj);
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
 * \brief    Compute internal concentrations
 * \details  Formula: rho*M*f
 * \param    void
 * \return   \e void
 */
void Model::compute_c( void )
{
  gsl_blas_dgemv(CblasNoTrans, _current_rho, _M, _f, 0.0, _c);
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
  _KM_f_product = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    _KM_f_product *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
  }
  gsl_vector_set(_tau_j, j, _KM_f_product/gsl_vector_get(_kcat_f, j));
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + inhibition (only one inhibitor per reaction)
 * \details  Formula: tau_j = prod(1+xc*1/KI[,j])*prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMMi( int j )
{
  _KM_f_product = 1.0;
  _KI_product   = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    _KM_f_product *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    _KI_product   *= 1.0+gsl_vector_get(_xc, i)*1.0/gsl_matrix_get(_KI, i, j);
  }
  gsl_vector_set(_tau_j, j, _KM_f_product*_KI_product/gsl_vector_get(_kcat_f, j));
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + activation (only one activator per reaction)
 * \details  Formula: tau_j = prod(1+KA[,j]/xc)*prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMMa( int j )
{
  _KM_f_product = 1.0;
  _KA_product   = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    _KM_f_product *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    _KA_product   *= 1.0+gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
  }
  gsl_vector_set(_tau_j, j, _KM_f_product*_KA_product/gsl_vector_get(_kcat_f, j));
}

/**
 * \brief    Irreversible Michaelis-Menten kinetics + inhibition + activation
 * \details  Formula: tau_j = prod(1+xc*1/KI[,j])*prod(1+KA[,j]/xc)*prod(1+Km_f[,j]/xc)/kcat_f[j]
 * \param    int j
 * \return   \e void
 */
void Model::iMMia( int j )
{
  _KM_f_product = 1.0;
  _KI_product   = 1.0;
  _KA_product   = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    _KM_f_product *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    _KI_product   *= 1.0+gsl_vector_get(_xc, i)*1.0/gsl_matrix_get(_KI, i, j);
    _KA_product   *= 1.0+gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
  }
  gsl_vector_set(_tau_j, j, _KM_f_product*_KI_product*_KA_product/gsl_vector_get(_kcat_f, j));
}

/**
 * \brief    Reversible Michaelis-Menten kinetics
 * \details  Formula: tau_j = 1/[ kcat_f[j]/prod(1+Km_f[,j]/xc) - kcat_b[j]/prod(1+Km_b[,j]/xc) ]
 * \param    int j
 * \return   \e void
 */
void Model::rMM( int j )
{
  _KM_f_product = 1.0;
  _KM_b_product = 1.0;
  for (int i = 0; i < _ni; i++)
  {
    _KM_f_product *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    _KM_b_product *= 1.0+gsl_matrix_get(_KM_b, i, j)/gsl_vector_get(_xc, i);
  }
  _KM_f_product = gsl_vector_get(_kcat_f, j)/_KM_f_product;
  _KM_b_product = gsl_vector_get(_kcat_b, j)/_KM_b_product;
  gsl_vector_set(_tau_j, j, 1/(_KM_f_product-_KM_b_product));
}

/**
 * \brief    Genome-scale Michaelis-Menten kinetics
 * \details  Irreversible MM considering reaction directionality
 * \param    int j
 * \return   \e void
 */
void Model::GMMM( int j )
{
  if (_directions[j] == "forward")
  {
    _KM_f_product = 1.0;
    for (int i = 0; i < _ni; i++)
    {
      _KM_f_product *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    }
    gsl_vector_set(_tau_j, j, _KM_f_product/gsl_vector_get(_kcat_f, j));
  }
  else if (_directions[j] == "backward")
  {
    _KM_b_product = 1.0;
    for (int i = 0; i < _ni; i++)
    {
      _KM_b_product *= 1.0+gsl_matrix_get(_KM_b, i, j)/gsl_vector_get(_xc, i);
    }
    gsl_vector_set(_tau_j, j, -_KM_b_product/gsl_vector_get(_kcat_b, j));
  }
  else if (_directions[j] == "reversible")
  {
    rMM(j);
  }
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
  /*-----------------------------*/
  /* 1) Define constants         */
  /*-----------------------------*/
  double term3 = gsl_vector_get(_kcat_f, j);
  /*-----------------------------*/
  /* 2) Calculate the derivative */
  /*-----------------------------*/
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
      double detivative = -term1*term2/term3;
      gsl_matrix_set(_ditau_j, j, i, detivative);
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
  /*-----------------------------*/
  /* 1) Define constants         */
  /*-----------------------------*/
  double constant_1 = 1.0;
  double constant_2 = 1.0;
  double term5      = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    constant_1 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    constant_2 *= 1.0+gsl_vector_get(_xc, i)*1.0/gsl_matrix_get(_KI, i, j);
  }
  /*-----------------------------*/
  /* 2) Calculate the derivative */
  /*-----------------------------*/
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double term1 = 1.0/gsl_matrix_get(_KI, y, j)*constant_1;
    double term2 = constant_2;
    double term3 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term4 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term4 *= 1.0+gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
      double derivative = term1-term2*term3*term4/term5;
      gsl_matrix_set(_ditau_j, j, i, derivative);
    }
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
  /*-----------------------------*/
  /* 1) Define constants         */
  /*-----------------------------*/
  double constant_1 = 1.0;
  double constant_2 = 1.0;
  double term6      = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    constant_1 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    constant_2 *= 1.0+gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
  }
  /*-----------------------------*/
  /* 2) Calculate the derivative */
  /*-----------------------------*/
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double term1 = gsl_matrix_get(_KA, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term2 = constant_1;
    double term3 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term4 = constant_2;
    double term5 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term5 *= 1.0+gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
      double derivative = -(term1*term2+term3*term4*term5)/term6;
      gsl_matrix_set(_ditau_j, j, i, derivative);
    }
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
  /*-----------------------------*/
  /* 1) Define constants         */
  /*-----------------------------*/
  double constant_1 = 1.0;
  double constant_2 = 1.0;
  double constant_3 = 1.0;
  double term9      = gsl_vector_get(_kcat_f, j);
  for (int i = 0; i < _ni; i++)
  {
    constant_1 *= 1.0+gsl_matrix_get(_KM_f, i, j)/gsl_vector_get(_xc, i);
    constant_2 *= 1.0+gsl_vector_get(_xc, i)*1.0/gsl_matrix_get(_KI, i, j);
    constant_3 *= 1.0+gsl_matrix_get(_KA, i, j)/gsl_vector_get(_xc, i);
  }
  /*-----------------------------*/
  /* 2) Calculate the derivative */
  /*-----------------------------*/
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double term1 = 1.0/gsl_matrix_get(_KI, y, j)*constant_1*constant_3;
    double term2 = constant_2;
    double term3 = -gsl_matrix_get(_KA, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term4 = constant_1;
    double term5 = constant_2;
    double term6 = constant_3;
    double term7 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
    double term8 = 1.0;
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term8 *= 1.0+gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
      }
      double derivative = term1+(term2*term3*term4)+(term5*term6*term7*term8)/term9;
      gsl_matrix_set(_ditau_j, j, i, derivative);
    }
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
  /*-----------------------------*/
  /* 1) Define constants         */
  /*-----------------------------*/
  double constant_1 = gsl_vector_get(_kcat_f, j);
  double constant_2 = gsl_vector_get(_kcat_b, j);
  /*-----------------------------*/
  /* 2) Calculate the derivative */
  /*-----------------------------*/
  for (int i = 0; i < _nc; i++)
  {
    int    y     = i+_nx;
    double term1 = 1.0;
    double term2 = gsl_matrix_get(_KM_f, y, j)/gsl_pow_int(gsl_matrix_get(_KM_f, y, j)+gsl_vector_get(_c, i), 2);
    double term3 = 1.0;
    double term4 = gsl_matrix_get(_KM_b, y, j)/gsl_pow_int(gsl_matrix_get(_KM_b, y, j)+gsl_vector_get(_c, i), 2);
    for (int index = 0; index < _ni; index++)
    {
      if (index != y)
      {
        term1 *= 1.0+gsl_matrix_get(_KM_f, index, j)/gsl_vector_get(_xc, index);
        term3 *= 1.0+gsl_matrix_get(_KM_b, index, j)/gsl_vector_get(_xc, index);
      }
      term1             = constant_1/term1;
      term3             = constant_2/term3;
      double detivative = (term1*term2-term3*term4)*(-gsl_vector_get(_tau_j, j));
      gsl_matrix_set(_ditau_j, j, i, detivative);
    }
  }
}

/**
 * \brief    Derivative of iMM with respect to metabolite concentrations, for genome-scale models
 * \details  Irreversible MM derivatives considering reaction directionality
 * \param    int j
 * \return   \e void
 */
void Model::dGMMM( int j )
{
  if (_directions[j] == "forward")
  {
    double term3 = gsl_vector_get(_kcat_f, j);
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
        double detivative = -term1*term2/term3;
        gsl_matrix_set(_ditau_j, j, i, detivative);
      }
    }
  }
  else if (_directions[j] == "backward")
  {
    double term3 = gsl_vector_get(_kcat_b, j);
    for (int i = 0; i < _nc; i++)
    {
      int    y     = i+_nx;
      double term1 = gsl_matrix_get(_KM_b, y, j)/gsl_pow_int(gsl_vector_get(_c, i), 2);
      double term2 = 1.0;
      for (int index = 0; index < _ni; index++)
      {
        if (index != y)
        {
          term2 *= 1.0+gsl_matrix_get(_KM_b, index, j)/gsl_vector_get(_xc, index);
        }
        double detivative = term1*term2/term3;
        gsl_matrix_set(_ditau_j, j, i, detivative);
      }
    }
  }
  else if (_directions[j] == "reversible")
  {
    drMM(j);
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
  gsl_blas_ddot(_tau_j, _f, &_mu);
  _mu = gsl_matrix_get(_M, _a, _r)*gsl_vector_get(_f, _r)/_mu;
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
  /*--------------------------------------------*/
  /* 1) Test density constraint                 */
  /*--------------------------------------------*/
  bool test1 = true;
  if (fabs(_density-1.0) >= DENSITY_CONSTRAINT_TOL)
  {
    test1 = false;
  }
  /*--------------------------------------------*/
  /* 2) Test negative concentrations constraint */
  /*--------------------------------------------*/
  bool test2 = true;
  if (gsl_vector_min(_c) < -NEGATIVE_C_TOL)
  {
    test2 = false;
  }
  /*--------------------------------------------*/
  /* 3) Test negative proteins constraint       */
  /*--------------------------------------------*/
  bool test3 = true;
  if (gsl_vector_min(_p) < -NEGATIVE_P_TOL)
  {
    test3 = false;
  }
  /*--------------------------------------------*/
  /* 4) Print error message if inconsistent     */
  /*--------------------------------------------*/
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
    /*** CASE 1: Reaction is irreversible and positive ***/
    if (_directions[j+1] == "forward" && gsl_vector_get(_f_trunc, j) < MIN_FLUX_FRACTION && gsl_vector_get(_GCC_f, j+1) < 0.0)
    {
      gsl_vector_set(_GCC_f, j+1, 0.0);
      gsl_vector_set(_f_trunc, j, MIN_FLUX_FRACTION);
    }
    /*** CASE 2: Reaction is irreversible and negative ***/
    else if (_directions[j+1] == "backward" && gsl_vector_get(_f_trunc, j) > -MIN_FLUX_FRACTION && gsl_vector_get(_GCC_f, j+1) > 0.0)
    {
      gsl_vector_set(_GCC_f, j+1, 0.0);
      gsl_vector_set(_f_trunc, j, -MIN_FLUX_FRACTION);
    }
    /*** CASE 3: Reaction is reversible and tends to zero ***/
    else if (_directions[j+1] == "reversible" && fabs(gsl_vector_get(_f_trunc, j)) < MIN_FLUX_FRACTION)
    {
      gsl_vector_set(_GCC_f, j+1, 0.0);
      if (gsl_vector_get(_f_trunc, j) > 0.0)
      {
        gsl_vector_set(_f_trunc, j, MIN_FLUX_FRACTION);
      }
      else if (gsl_vector_get(_f_trunc, j) < 0.0)
      {
        gsl_vector_set(_f_trunc, j, -MIN_FLUX_FRACTION);
      }
    }
  }
}

/**
 * \brief    Open trajectory output files
 * \details  Also writes headers
 * \param    void
 * \return   \e void
 */
void Model::open_trajectory_output_files( void )
{
  /*------------------*/
  /* 1) Open files    */
  /*------------------*/
  std::stringstream model_trajectory_filename;
  std::stringstream f_trajectory_filename;
  std::stringstream c_trajectory_filename;
  std::stringstream v_trajectory_filename;
  std::stringstream p_trajectory_filename;
  std::stringstream b_trajectory_filename;
  model_trajectory_filename << "./output/" << _model_name << "_model_trajectory.csv";
  f_trajectory_filename << "./output/" << _model_name << "_f_trajectory.csv";
  c_trajectory_filename << "./output/" << _model_name << "_c_trajectory.csv";
  v_trajectory_filename << "./output/" << _model_name << "_v_trajectory.csv";
  p_trajectory_filename << "./output/" << _model_name << "_p_trajectory.csv";
  b_trajectory_filename << "./output/" << _model_name << "_b_trajectory.csv";
  _model_trajectory_file.open(model_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _f_trajectory_file.open(f_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _c_trajectory_file.open(c_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _v_trajectory_file.open(v_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _p_trajectory_file.open(p_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  _b_trajectory_file.open(b_trajectory_filename.str(), std::ios::out | std::ios::trunc);
  /*------------------*/
  /* 2) Write headers */
  /*------------------*/
   _model_trajectory_file << "t;dt;mu;density;consistent\n";
  _f_trajectory_file << "t;dt";
  _c_trajectory_file << "t;dt";
  _v_trajectory_file << "t;dt";
  _p_trajectory_file << "t;dt";
  _b_trajectory_file << "t;dt";
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
 * \param    double t
 * \param    double dt
 * \return   \e void
 */
void Model::write_trajectory_output_files( double t, double dt )
{
  _model_trajectory_file << t << ";" << dt << ";" << _mu << ";" << _density << ";" << _consistent << "\n";
  _f_trajectory_file << t << ";" << dt;
  _c_trajectory_file << t << ";" << dt;
  _v_trajectory_file << t << ";" << dt;
  _p_trajectory_file << t << ";" << dt;
  _b_trajectory_file << t << ";" << dt;
  for (int i = 0; i < _nc; i++)
  {
    _c_trajectory_file << ";" << gsl_vector_get(_c, i);
    _b_trajectory_file << ";" << gsl_vector_get(_b, i);
  }
  for (int j = 0; j < _nj; j++)
  {
    _f_trajectory_file << ";" << gsl_vector_get(_f, j);
    _v_trajectory_file << ";" << gsl_vector_get(_v, j);
    _p_trajectory_file << ";" << gsl_vector_get(_p, j);
  }
  _f_trajectory_file << "\n";
  _c_trajectory_file << "\n";
  _v_trajectory_file << "\n";
  _p_trajectory_file << "\n";
  _b_trajectory_file << "\n";
}

/**
 * \brief    Close trajectory output files
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::close_trajectory_ouput_files( void )
{
  _model_trajectory_file.close();
  _f_trajectory_file.close();
  _c_trajectory_file.close();
  _v_trajectory_file.close();
  _p_trajectory_file.close();
  _b_trajectory_file.close();
}

/**
 * \brief    Open optimum output files
 * \details  Also writes headers
 * \param    void
 * \return   \e void
 */
void Model::open_optimum_output_files( void )
{
  /*------------------*/
  /* 1) Open files    */
  /*------------------*/
  std::stringstream model_optimum_filename;
  std::stringstream f_optimum_filename;
  std::stringstream c_optimum_filename;
  std::stringstream v_optimum_filename;
  std::stringstream p_optimum_filename;
  std::stringstream b_optimum_filename;
  model_optimum_filename << "./output/" << _model_name << "_model_optimum.csv";
  f_optimum_filename << "./output/" << _model_name << "_f_optimum.csv";
  c_optimum_filename << "./output/" << _model_name << "_c_optimum.csv";
  v_optimum_filename << "./output/" << _model_name << "_v_optimum.csv";
  p_optimum_filename << "./output/" << _model_name << "_p_optimum.csv";
  b_optimum_filename << "./output/" << _model_name << "_b_optimum.csv";
  _model_optimum_file.open(model_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _f_optimum_file.open(f_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _c_optimum_file.open(c_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _v_optimum_file.open(v_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _p_optimum_file.open(p_optimum_filename.str(), std::ios::out | std::ios::trunc);
  _b_optimum_file.open(b_optimum_filename.str(), std::ios::out | std::ios::trunc);
  /*------------------*/
  /* 2) Write headers */
  /*------------------*/
   _model_optimum_file << "condition;mu;density;consistent;converged\n";
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
 * \return   \e void
 */
void Model::write_optimum_output_files( std::string condition, bool converged )
{
  _model_optimum_file << condition << ";" << _mu << ";" << _density << ";" << _consistent << ";" << converged << "\n";
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
}

/**
 * \brief    Close optimum output files
 * \details  --
 * \param    void
 * \return   \e void
 */
void Model::close_optimum_ouput_files( void )
{
  _model_optimum_file.close();
  _f_optimum_file.close();
  _c_optimum_file.close();
  _v_optimum_file.close();
  _p_optimum_file.close();
  _b_optimum_file.close();
}

/**
 * \brief    Compute the gradient ascent trajectory for small models exclusively
 * \details  Basic trajectory controls are implemented
 * \param    std::string condition
 * \param    double initial_dt
 * \param    double max_t
 * \param    bool save_trajectory
 * \return   \e bool
 */
bool Model::compute_gradient_ascent_trajectory_for_small_models( std::string condition, double initial_dt, double max_t, bool save_trajectory )
{
  assert(_condition_indices.find(condition) != _condition_indices.end());
  assert(initial_dt > 0.0);
  assert(max_t > 0.0);
  if (save_trajectory)
  {
    open_trajectory_output_files();
  }
  set_condition(condition);
  initialize_f();
  calculate();
  if (!_consistent)
  {
    throw std::invalid_argument("> Model initial state f0 is inconsistent");
  }
  gsl_vector* previous_f_trunc = gsl_vector_alloc(_nj-1);
  gsl_vector* scaled_dmudt     = gsl_vector_alloc(_nj-1);
  gsl_vector_memcpy(previous_f_trunc, _f_trunc);
  double previous_mu         = 0.0;
  double t                   = 0.0;
  double dt                  = initial_dt;
  int    dt_counter          = 0;
  int    constant_mu_counter = 0;
  int    nb_iterations       = 0;
  int    nb_successes        = 0;
  while (t < max_t)
  {
    assert(dt > 1e-100);
    nb_iterations++;
    if (constant_mu_counter >= TRAJECTORY_STABLE_MU_COUNT)
    {
      break;
    }
    previous_mu           = _mu;
    gsl_vector_view dmudt = gsl_vector_subvector(_GCC_f, 1, _nj-1);
    gsl_vector_memcpy(scaled_dmudt, &dmudt.vector);
    gsl_vector_scale(scaled_dmudt, dt);
    gsl_vector_add(_f_trunc, scaled_dmudt);
    set_f();
    calculate();
    if (_consistent)
    {
      gsl_vector_memcpy(previous_f_trunc, _f_trunc);
      nb_successes++;
      t = t+dt;
      dt_counter++;
      if (save_trajectory)
      {
        write_trajectory_output_files(t, dt);
      }
      if (fabs(_mu-previous_mu) < TRAJECTORY_CONVERGENCE_TOL)
      {
        constant_mu_counter++;
      }
      else
      {
        constant_mu_counter = 0;
      }
      if (dt_counter == 1000)
      {
        dt         *= INCREASING_DT_FACTOR;
        dt_counter  = 0;
      }
    }
    else
    {
      gsl_vector_memcpy(_f_trunc, previous_f_trunc);
      set_f();
      calculate();
      assert(_consistent);
      dt         /= DECREASING_DT_FACTOR;
      dt_counter  = 0;
    }
  }
  gsl_vector_free(previous_f_trunc);
  gsl_vector_free(scaled_dmudt);
  previous_f_trunc = NULL;
  scaled_dmudt     = NULL;
  if (save_trajectory)
  {
    close_trajectory_ouput_files();
  }
  if (constant_mu_counter < TRAJECTORY_STABLE_MU_COUNT)
  {
    std::cout << "> Convergence not reached after T=" << max_t << " (nb iterations=" << nb_iterations << ")" << std::endl;
    return(false);
  }
  else
  {
    std::cout << "> Convergence reached (mu=" << _mu << ", nb iterations=" << nb_iterations << ")" << std::endl;
    return(true);
  }
}

/**
 * \brief    Compute the gradient ascent trajectory for genome-scale models
 * \details  More elaborated controls are implemented on the trajectory, and specific functions are used
 * \param    std::string condition
 * \param    double initial_dt
 * \param    double max_t
 * \param    bool save_trajectory
 * \return   \e bool
 */
bool Model::compute_gradient_ascent_trajectory_for_genome_scale_models( std::string condition, double initial_dt, double max_t, bool save_trajectory )
{
  assert(_condition_indices.find(condition) != _condition_indices.end());
  assert(initial_dt > 0.0);
  assert(max_t > 0.0);
  if (save_trajectory)
  {
    open_trajectory_output_files();
  }
  set_condition(condition);
  initialize_f();
  calculate_GM();
  if (!_consistent)
  {
    throw std::invalid_argument("> Model initial state f0 is inconsistent");
  }
  gsl_vector* previous_f_trunc = gsl_vector_alloc(_nj-1);
  gsl_vector* scaled_dmudt     = gsl_vector_alloc(_nj-1);
  gsl_vector_memcpy(previous_f_trunc, _f_trunc);
  double previous_mu         = 0.0;
  double t                   = 0.0;
  double dt                  = initial_dt;
  int    dt_counter          = 0;
  int    constant_mu_counter = 0;
  int    nb_iterations       = 0;
  int    nb_successes        = 0;
  while (t < max_t)
  {
    /*------------------------------------------------*/
    /* 1) Check the size of dt                        */
    /*------------------------------------------------*/
    if (dt <= 1e-100)
    {
      throw std::invalid_argument("> dt is too small");
    }
    nb_iterations++;
    /*------------------------------------------------*/
    /* 2) Check trajectory convergence                */
    /*------------------------------------------------*/
    if (constant_mu_counter >= TRAJECTORY_STABLE_MU_COUNT)
    {
      break;
    }
    set_f();
    calculate_second_order_GM();
    block_reactions();
    previous_mu           = _mu;
    gsl_vector_view dmudt = gsl_vector_subvector(_GCC_f, 1, _nj-1);
    gsl_vector_memcpy(scaled_dmudt, &dmudt.vector);
    gsl_vector_scale(scaled_dmudt, dt);
    gsl_vector_add(_f_trunc, scaled_dmudt);
    set_f();
    calculate_first_order_GM();
    if (_consistent)
    {
      gsl_vector_memcpy(previous_f_trunc, _f_trunc);
      nb_successes++;
      t = t+dt;
      dt_counter++;
      if (save_trajectory)
      {
        write_trajectory_output_files(t, dt);
      }
      if (fabs(_mu-previous_mu) < TRAJECTORY_CONVERGENCE_TOL)
      {
        constant_mu_counter++;
      }
      else
      {
        constant_mu_counter = 0;
      }
      if (dt_counter == 1000)
      {
        dt         *= INCREASING_DT_FACTOR;
        dt_counter  = 0;
      }
    }
    else
    {
      gsl_vector_memcpy(_f_trunc, previous_f_trunc);
      set_f();
      calculate_GM();
      assert(_consistent);
      dt         /= DECREASING_DT_FACTOR;
      dt_counter  = 0;
    }
  }
  gsl_vector_free(previous_f_trunc);
  gsl_vector_free(scaled_dmudt);
  previous_f_trunc = NULL;
  scaled_dmudt     = NULL;
  if (save_trajectory)
  {
    close_trajectory_ouput_files();
  }
  if (constant_mu_counter < TRAJECTORY_STABLE_MU_COUNT)
  {
    std::cout << "> Convergence not reached after T=" << max_t << " (nb iterations=" << nb_iterations << ")" << std::endl;
    return(false);
  }
  else
  {
    std::cout << "> Convergence reached (mu=" << _mu << ", nb iterations=" << nb_iterations << ")" << std::endl;
    return(true);
  }
}


