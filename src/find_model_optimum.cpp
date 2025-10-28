/**
 * \file      find_model_optimum.cpp
 * \author    Charles Rocabert
 * \date      22-07-2024
 * \copyright gbacpp. Copyright © 2024-2025 Charles Rocabert.
 * \license   GNU General Public License v3 (GPLv3)
 * \brief     find_model_optimum executable
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <ctime>
#include <vector>
#include <unordered_map>
#include <sys/stat.h>
#include <assert.h>

#include "Config.h"
//#include "../cmake/Config.h"
#include "./lib/Macros.hpp"
#include "./lib/Enums.hpp"
#include "./lib/Model.hpp"

void readArgs( int argc, char const** argv, std::string &model_path, std::string &model_name, std::string &condition, bool &print_optimum, bool &write_optimum, bool &write_trajectory, std::string &output_path, double &tol, double &mu_tol, double &q_tol, int &convergence_count, int &max_iter, bool &hessian, bool &reload, bool &restart, bool &use_previous_sol, bool &verbose, bool &extra_verbose );
void printUsage( void );
void printHeader( void );


/**
 * \brief    main function
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \return   \e int
 */
int main(int argc, char const** argv)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Read parameters                                */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  std::string model_path        = "";
  std::string model_name        = "";
  std::string condition         = "";
  bool        print_optimum     = false;
  bool        write_optimum     = false;
  bool        write_trajectory  = false;
  std::string output_path       = ".";
  double      tol               = 1e-10;
  double      mu_tol            = 1e-10;
  double      q_tol             = 1e-10;
  int         convergence_count = 10000;
  int         max_iter          = 100000000;
  bool        hessian           = false;
  bool        reload            = false;
  bool        restart           = false;
  bool        use_previous_sol  = false;
  bool        verbose           = false;
  bool        extra_verbose     = false;
  readArgs(argc, argv, model_path, model_name, condition, print_optimum, write_optimum, write_trajectory, output_path, tol, mu_tol, q_tol, convergence_count, max_iter, hessian, reload, restart, use_previous_sol, verbose, extra_verbose);
  if (condition != "all" && use_previous_sol)
  {
    throw std::invalid_argument("> Error: option -previous (--use-previous-sol) can only be used with condition \"all\"");
  }
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Print the header in verbose mode               */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  if ((verbose || extra_verbose) && condition != "all")
  {
    printHeader();
    std::cout << "OPTIMIZING MODEL \"" << model_name << "\" FOR CONDITION " << condition << ":" << std::endl;
  }
  if ((verbose || extra_verbose) && condition == "all")
  {
    printHeader();
    std::cout << "OPTIMIZING MODEL \"" << model_name << "\" FOR ALL CONDITIONS:" << std::endl;
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Load the model                                 */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  Model* model = new Model(model_path, model_name);
  model->read_from_csv();
  model->set_tol(tol);
  model->set_mu_tol(mu_tol);
  model->set_q_tol(q_tol);
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) Run the calculation depending on the condition */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  if (condition != "all")// && condition != "random")
  {
    model->compute_optimum(condition, print_optimum, write_optimum, write_trajectory, output_path, convergence_count, max_iter, hessian, reload, restart, verbose, extra_verbose);
  }
  else
  {
    model->compute_optimum_by_condition(print_optimum, write_optimum, write_trajectory, output_path, convergence_count, max_iter, hessian, reload, restart, use_previous_sol, verbose, extra_verbose);
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 5) Free memory and exit                           */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  delete model;
  model = NULL;
  return EXIT_SUCCESS;
}

/**
 * \brief    Read command line arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    std::string &model_path
 * \param    std::string &model_name
 * \param    std::string &condition
 * \param    bool &print_optimum
 * \param    bool &write_optimum
 * \param    bool &write_trajectory
 * \param    std::string &output_path
 * \param    double &tol
 * \param    double &mu_tol
 * \param    double &q_tol
 * \param    double &convergence_count
 * \param    int &max_iter
 * \param    bool &hessian
 * \param    bool &reload
 * \param    bool &restart
 * \param    bool &use_previous_sol
 * \param    bool &verbose
 * \param    bool &extra_verbose
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string &model_path, std::string &model_name, std::string &condition, bool &print_optimum, bool &write_optimum, bool &write_trajectory, std::string &output_path, double &tol, double &mu_tol, double &q_tol, int &convergence_count, int &max_iter, bool &hessian, bool &reload, bool &restart, bool &use_previous_sol, bool &verbose, bool &extra_verbose )
{
  if (argc == 1)
  {
    printf("You must provide all the mandatory arguments (see -h or --help). Exit.\n");
    exit(EXIT_SUCCESS);
  }
  int counter = 0;
  for (int i = 0; i < argc; i++)
  {
    /****************************************************************/
    
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
    {
      printUsage();
      exit(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-version") == 0 || strcmp(argv[i], "--version") == 0)
    {
      std::cout << PACKAGE << " (" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << ")\n";
      exit(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-path") == 0 || strcmp(argv[i], "--model-path") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: Model path value is missing");
      }
      else
      {
        model_path = argv[i+1];
        counter++;
      }
    }
    else if (strcmp(argv[i], "-name") == 0 || strcmp(argv[i], "--model-name") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: Model name value is missing");
      }
      else
      {
        model_name = argv[i+1];
        counter++;
      }
    }
    else if (strcmp(argv[i], "-condition") == 0 || strcmp(argv[i], "--condition") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: Condition is missing");
      }
      else
      {
        condition = argv[i+1];
        counter++;
      }
    }
    else if (strcmp(argv[i], "-print") == 0 || strcmp(argv[i], "--print-trajectory") == 0)
    {
      print_optimum = true;
    }
    else if (strcmp(argv[i], "-optimum") == 0 || strcmp(argv[i], "--write-optimum") == 0)
    {
      write_optimum = true;
    }
    else if (strcmp(argv[i], "-trajectory") == 0 || strcmp(argv[i], "--write-trajectory") == 0)
    {
      write_trajectory = true;
    }
    else if (strcmp(argv[i], "-output") == 0 || strcmp(argv[i], "--output-path") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: Output path value is missing");
      }
      else
      {
        output_path = argv[i+1];
      }
    }
    else if (strcmp(argv[i], "-tol") == 0 || strcmp(argv[i], "--tolerance") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: Tolerance value is missing");
      }
      else
      {
        tol = atof(argv[i+1]);
      }
    }
    else if (strcmp(argv[i], "-mutol") == 0 || strcmp(argv[i], "--mu-tolerance") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: Mu tolerance value is missing");
      }
      else
      {
        mu_tol = atof(argv[i+1]);
      }
    }
    else if (strcmp(argv[i], "-qtol") == 0 || strcmp(argv[i], "--q-tolerance") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: q tolerance value is missing");
      }
      else
      {
        q_tol = atof(argv[i+1]);
      }
    }
    else if (strcmp(argv[i], "-conv") == 0 || strcmp(argv[i], "--convergence-count") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: Mu convergence count is missing");
      }
      else
      {
        convergence_count = atof(argv[i+1]);
      }
    }
    else if (strcmp(argv[i], "-max") == 0 || strcmp(argv[i], "--max-iter") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: maximum number of iterations is missing");
      }
      else
      {
        max_iter = atoi(argv[i+1]);
      }
    }
    else if (strcmp(argv[i], "-hessian") == 0 || strcmp(argv[i], "--hessian") == 0)
    {
      hessian = true;
    }
    else if (strcmp(argv[i], "-reload") == 0 || strcmp(argv[i], "--reload") == 0)
    {
      reload = true;
    }
    else if (strcmp(argv[i], "-restart") == 0 || strcmp(argv[i], "--restart") == 0)
    {
      restart = true;
    }
    else if (strcmp(argv[i], "-previous") == 0 || strcmp(argv[i], "--use-previous-sol") == 0)
    {
      use_previous_sol = true;
    }
    else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0)
    {
      verbose = true;
    }
    else if (strcmp(argv[i], "-vv") == 0 || strcmp(argv[i], "--extra-verbose") == 0)
    {
      extra_verbose = true;
    }
  }
  if (counter < 3)
  {
    throw std::invalid_argument("> You must provide all the mandatory arguments (see -h or --help)");
  }
}

/**
 * \brief    Print usage
 * \details  --
 * \param    void
 * \return   \e void
 */
void printUsage( void )
{
  std::cout << "\n";
  std::cout << "************************************************************************\n";
#ifdef DEBUG
  std::cout << "* " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " (Debug)\n";
#endif
#ifdef NDEBUG
  std::cout << "* " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " (Release)\n";
#endif
  std::cout << "* Web: https://github.com/charlesrocabert/gbacpp\n";
  std::cout << "* Copyright © 2024-2025 Charles Rocabert.\n";
  std::cout << "*\n";
  std::cout << "* This program is free software: you can redistribute it and/or modify\n";
  std::cout << "* it under the terms of the GNU General Public License as published by\n";
  std::cout << "* the Free Software Foundation, either version 3 of the License, or\n";
  std::cout << "* (at your option) any later version.\n";
  std::cout << "*\n";
  std::cout << "* This program is distributed in the hope that it will be useful,\n";
  std::cout << "* but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
  std::cout << "* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
  std::cout << "* GNU General Public License for more details.\n";
  std::cout << "*\n";
  std::cout << "* You should have received a copy of the GNU General Public License\n";
  std::cout << "* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n";
  std::cout << "************************************************************************\n";
  std::cout << "\n";
  std::cout << "Usage: find_model_optimum -h or --help\n";
  std::cout << "   or: find_model_optimum [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -version, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -path, --model-path (MANDATORY)\n";
  std::cout << "        specify the path of the model to be loaded\n";
  std::cout << "  -name, --model-name (MANDATORY)\n";
  std::cout << "        specify the name of the model to be loaded\n";
  std::cout << "  -condition, --condition (MANDATORY)\n";
  std::cout << "        specify the condition (condition identifier / all)\n";
  std::cout << "  -print, --print-optimum\n";
  std::cout << "        indicates if the optimum should be printed in the standard output\n";
  std::cout << "  -optimum, --write-optimum\n";
  std::cout << "        indicates if the optimum should be written in output files\n";
  std::cout << "  -trajectory, --write-trajectory\n";
  std::cout << "        indicates if the trajectory should be written in output files\n";
  std::cout << "  -output, --output-path\n";
  std::cout << "        specify the path of output files\n";
  std::cout << "  -tol, --tolerance\n";
  std::cout << "        specify the tolerance value\n";
  std::cout << "  -mutol, --mu-tolerance\n";
  std::cout << "        specify the relative growth rate difference tolerance value to assume convergence\n";
  std::cout << "  -qtol, --q-tolerance\n";
  std::cout << "        specify the maximal relative q difference tolerance value to assume convergence\n";
  std::cout << "  -conv, --convergence-count\n";
  std::cout << "        specify the number of iterations under mu tolerance needed to assume convergence\n";
  std::cout << "  -max, --max-iter\n";
  std::cout << "        specify the maximal number of iterations\n";
  //std::cout << "  -hessian, --hessian\n";
  //std::cout << "        indicates if the diagonal Hessian should be estimated\n";
  //std::cout << "  -reload, --reload\n";
  //std::cout << "        indicates if the last trajectory point should be used as q0\n";
  //std::cout << "  -restart, --restart\n";
  //std::cout << "        indicates if the last trajectory point should be used as a fresh start\n";
  std::cout << "  -previous, --use-previous-sol\n";
  std::cout << "        indicates if the last optimal solution should be used (only for condition=\"all\")\n";
  std::cout << "  -v, --verbose\n";
  std::cout << "        indicates if the program should run in verbose mode\n";
  std::cout << "  -vv, --extra-verbose\n";
  std::cout << "        indicates if the program should run in extra-verbose mode\n";
  std::cout << "\n";
}

/**
 * \brief    Print header
 * \details  --
 * \param    void
 * \return   \e void
 */
void printHeader( void )
{
  std::cout << "\n";
  std::cout << "************************************************************************\n";
#ifdef DEBUG
  std::cout << "* " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " (Debug)\n";
#endif
#ifdef NDEBUG
  std::cout << "* " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " (Release)\n";
#endif
  std::cout << "* Web: https://github.com/charlesrocabert/gbacpp                        \n";
  std::cout << "* GPLv3 License © 2024-2025 Charles Rocabert.                           \n";
  std::cout << "************************************************************************\n";
  std::cout << "\n";
}

