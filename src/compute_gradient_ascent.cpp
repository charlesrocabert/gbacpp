/**
 * \file      compute_gradient_ascent.cpp
 * \author    Charles Rocabert
 * \date      22-07-2024
 * \copyright GBAcpp. Copyright © 2024-2025 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     compute_gradient_ascent executable
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

#include "../cmake/Config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <ctime>
#include <vector>
#include <unordered_map>
#include <sys/stat.h>
#include <assert.h>

#include "./lib/Macros.hpp"
#include "./lib/Enums.hpp"
#include "./lib/Model.hpp"

void readArgs( int argc, char const** argv, std::string &path, std::string &name, std::string &condition, double &initial_dt, double &max_t, int &max_mu_count, bool &save, std::string &output_path );
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
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Read parameters                             */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  std::string path         = "";
  std::string name         = "";
  std::string condition    = "";
  double      initial_dt   = 0.0;
  double      max_t        = 0.0;
  int         max_mu_count = 0;
  bool        save         = false;
  std::string output       = "";
  readArgs(argc, argv, path, name, condition, initial_dt, max_t, max_mu_count, save, output);
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Load the model and calculate the trajectory */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  std::clock_t begin = clock();
  Model*       model     = new Model(path, name);
  model->set_condition(condition);
  bool         converged = model->compute_gradient_ascent(condition, initial_dt, max_t, max_mu_count, save, output);
  std::clock_t end       = clock();
  double       runtime   = double(end-begin)/CLOCKS_PER_SEC;
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Save the optimum in case of convergence     */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  model->open_optimum_output_files(output, condition);
  model->write_optimum_output_files(condition, converged, runtime);
  model->close_optimum_ouput_files();
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) Free memory and exit                        */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  delete model;
  model = NULL;
  std::cout << "Elapsed time: " << runtime << " seconds" << std::endl;
  return EXIT_SUCCESS;
}

/**
 * \brief    Read command line arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    std::string &path
 * \param    std::string &name
 * \param    std::string &condition
 * \param    double &initial_dt
 * \param    double &max_t
 * \param    int &max_mu_count
 * \param    bool &save
 * \param    std::string &output
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string &path, std::string &name, std::string &condition, double &initial_dt, double &max_t, int &max_mu_count, bool &save, std::string &output )
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
        throw std::invalid_argument("> model path value is missing");
      }
      else
      {
        path = argv[i+1];
        counter++;
      }
    }
    else if (strcmp(argv[i], "-name") == 0 || strcmp(argv[i], "--model-name") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> model name value is missing");
      }
      else
      {
        name = argv[i+1];
        counter++;
      }
    }
    else if (strcmp(argv[i], "-condition") == 0 || strcmp(argv[i], "--condition") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> condition value is missing");
      }
      else
      {
        condition = argv[i+1];
        counter++;
      }
    }
    else if (strcmp(argv[i], "-dt") == 0 || strcmp(argv[i], "--initial-dt") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> initial dt value is missing");
      }
      else
      {
        initial_dt = atof(argv[i+1]);
        counter++;
      }
    }
    else if (strcmp(argv[i], "-maxt") == 0 || strcmp(argv[i], "--max-time") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> max time value is missing");
      }
      else
      {
        max_t = atof(argv[i+1]);
        counter++;
      }
    }
    else if (strcmp(argv[i], "-max-mu-count") == 0 || strcmp(argv[i], "--max-mu-count") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> max mu count value is missing");
      }
      else
      {
        max_mu_count = atoi(argv[i+1]);
      }
    }
    else if (strcmp(argv[i], "-save") == 0 || strcmp(argv[i], "--save-trajectory") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> output path value is missing");
      }
      else
      {
        save = true;
      }
    }
    else if (strcmp(argv[i], "-output") == 0 || strcmp(argv[i], "--output-path") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> output path value is missing");
      }
      else
      {
        output = argv[i+1];
      }
    }
  }
  if (counter < 5)
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
  std::cout << " GBAcpp (Growth Balance Analysis for C++)                            \n";
  std::cout << " Copyright © 2024-2025 Charles Rocabert                              \n";
  std::cout << " Web: https://github.com/charlesrocabert/gbacpp                      \n";
  std::cout << "                                                                     \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                     \n";
  std::cout << " This is free software, and you are welcome to redistribute it under \n";
  std::cout << " certain conditions; See the GNU General Public License for details  \n";
  std::cout << "*********************************************************************\n";
  std::cout << "\n";
  std::cout << "Usage: compute_gradient_ascent -h or --help\n";
  std::cout << "   or: compute_gradient_ascent [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -path, --model-path (MANDATORY)\n";
  std::cout << "        specify the path of the GBA model to be loaded\n";
  std::cout << "  -name, --model-name (MANDATORY)\n";
  std::cout << "        specify the name of the GBA model to be loaded\n";
  std::cout << "  -condition, --condition (MANDATORY)\n";
  std::cout << "        specify the external condition identifier\n";
  std::cout << "  -dt, --initial-dt (MANDATORY)\n";
  std::cout << "        specify the initial gradient timestep\n";
  std::cout << "  -maxt, --max-time (MANDATORY)\n";
  std::cout << "        specify the maximal time\n";
  std::cout << "  -max-mu-count, --max-mu-count\n";
  std::cout << "        specify the maximal number of iterations with unchanged mu\n";
  std::cout << "  -save, --save-trajectory\n";
  std::cout << "        specify if the trajectory should be saved as output files\n";
  std::cout << "  -output, --output-path\n";
  std::cout << "        specify the path of output files\n";
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
  std::cout << " GBAcpp (Growth Balance Analysis for C++)                            \n";
  std::cout << " Copyright © 2024-2025 Charles Rocabert                              \n";
  std::cout << " Web: https://github.com/charlesrocabert/gbacpp                      \n";
  std::cout << "                                                                     \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                     \n";
  std::cout << " This is free software, and you are welcome to redistribute it under \n";
  std::cout << " certain conditions; See the GNU General Public License for details  \n";
  std::cout << "*********************************************************************\n";
  std::cout << "\n";
}


