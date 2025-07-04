/**
 * \file      compute_optimum.cpp
 * \author    Charles Rocabert
 * \date      22-07-2024
 * \copyright gbacpp. Copyright © 2024-2025 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     compute_optimum executable
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
#include "./lib/Macros.hpp"
#include "./lib/Enums.hpp"
#include "./lib/Model.hpp"

void readArgs( int argc, char const** argv, std::string &model_path, std::string &model_name, std::string &condition, bool &print_optimum, bool &write_trajectory, std::string &output_path, double &tol, int &stable_count, double &max_time, bool &verbose );
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
  std::string model_path       = "";
  std::string model_name       = "";
  std::string condition        = "";
  bool        print_optimum    = false;
  bool        write_trajectory = false;
  std::string output_path      = "";
  double      tol              = 1e-10;
  int         stable_count     = 10000;
  double      max_time         = 100000.0;
  bool        verbose          = false;
  readArgs(argc, argv, model_path, model_name, condition, print_optimum, write_trajectory, output_path, tol, stable_count, max_time, verbose);
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Load the model                                 */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  Model* model = new Model(model_path, model_name);
  model->set_tol(tol);
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Run the calculation depending on the condition */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  if (condition != "all")// && condition != "random")
  {
    model->compute_optimum(condition, print_optimum, write_trajectory, output_path, stable_count, max_time, verbose);
  }
  else
  {
    model->compute_optimum_by_condition(print_optimum, write_trajectory, output_path, stable_count, max_time, verbose);
  }
  /*
  else if (condition == "random")
  {
    model->read_random_solutions();
    model->compute_optimum_by_random_solution(condition, print_optimum, write_trajectory, output_path, stable_count, max_time, verbose);
  }
  */
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) Free memory and exit                           */
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
 * \param    bool &write_trajectory
 * \param    std::string &output_path
 * \param    double &tol
 * \param    double &stable_count
 * \param    double &max_time
 * \param    bool &verbose
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string &model_path, std::string &model_name, std::string &condition, bool &print_optimum, bool &write_trajectory, std::string &output_path, double &tol, int &stable_count, double &max_time, bool &verbose )
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
    else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0)
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
    else if (strcmp(argv[i], "-write") == 0 || strcmp(argv[i], "--write-trajectory") == 0)
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
    else if (strcmp(argv[i], "-stable") == 0 || strcmp(argv[i], "--stable-count") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: Stable mu count is missing");
      }
      else
      {
        stable_count = atof(argv[i+1]);
      }
    }
    else if (strcmp(argv[i], "-max") == 0 || strcmp(argv[i], "--max-time") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> Error: Trajectory max time is missing");
      }
      else
      {
        max_time = atof(argv[i+1]);
      }
    }
    else if (strcmp(argv[i], "-verbose") == 0 || strcmp(argv[i], "--verbose") == 0)
    {
      verbose = true;
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
  std::cout << "*********************************************************************\n";
#ifdef DEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( debug )\n";
#endif
#ifdef NDEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( release )\n";
#endif
  std::cout << "                                                                     \n";
  std::cout << " gbacpp (growth balance analysis for C++)                            \n";
  std::cout << " Copyright © 2024-2025 Charles Rocabert                              \n";
  std::cout << " Web: https://github.com/charlesrocabert/gbacpp                      \n";
  std::cout << "                                                                     \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                     \n";
  std::cout << " This is free software, and you are welcome to redistribute it under \n";
  std::cout << " certain conditions; See the GNU General Public License for details  \n";
  std::cout << "*********************************************************************\n";
  std::cout << "\n";
  std::cout << "Usage: compute_optimum -h or --help\n";
  std::cout << "   or: compute_optimum [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -path, --model-path (MANDATORY)\n";
  std::cout << "        specify the path of the CGM to be loaded\n";
  std::cout << "  -name, --model-name (MANDATORY)\n";
  std::cout << "        specify the name of the CGM to be loaded\n";
  std::cout << "  -condition, --condition (MANDATORY)\n";
  std::cout << "        specify the condition (condition identifier / all)\n";
  std::cout << "  -print, --print-optimum\n";
  std::cout << "        indicates if the optimum should be printed in the standard output\n";
  std::cout << "  -write, --write-trajectory\n";
  std::cout << "        indicates if the trajectory should be written in output files\n";
  std::cout << "  -output, --output-path\n";
  std::cout << "        specify the path of output files\n";
  std::cout << "  -tol, --tolerance\n";
  std::cout << "        specify the tolerance value\n";
  std::cout << "  -stable-count, --stable-count\n";
  std::cout << "        specify the maximal number of iterations with unchanged mu\n";
  std::cout << "  -maxt, --max-time\n";
  std::cout << "        specify the maximal trajectory time\n";
  std::cout << "  -verbose, --verbose\n";
  std::cout << "        indicates if the program should run in verbose mode\n";
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
  std::cout << "*********************************************************************\n";
#ifdef DEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( debug )\n";
#endif
#ifdef NDEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( release )\n";
#endif
  std::cout << "                                                                     \n";
  std::cout << " gbacpp (growth balance analysis for C++)                           \n";
  std::cout << " Copyright © 2024-2025 Charles Rocabert                              \n";
  std::cout << " Web: https://github.com/charlesrocabert/gbacpp                      \n";
  std::cout << "                                                                     \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                     \n";
  std::cout << " This is free software, and you are welcome to redistribute it under \n";
  std::cout << " certain conditions; See the GNU General Public License for details  \n";
  std::cout << "*********************************************************************\n";
  std::cout << "\n";
}

