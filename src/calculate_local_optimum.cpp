/**
 * \file      calculate_local_optimum.cpp
 * \author    Charles Rocabert
 * \date      29-08-2024
 * \copyright GBA_Evolution. Copyright © 2024 Charles Rocabert. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     calculate_local_optimum executable
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

#include "../cmake/Config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <sys/stat.h>
#include <assert.h>

#include "./lib/Macros.hpp"
#include "./lib/Enums.hpp"
#include "./lib/Structs.hpp"
#include "./lib/Model.hpp"

void readArgs( int argc, char const** argv, std::string &path, std::string &name, msize &size, double &initial_dt, double &max_t, bool &save, std::string &output_path, bool &parallel_computing );
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
  std::clock_t begin = clock();
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Read parameters                             */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  std::string path       = "";
  std::string name       = "";
  msize       size       = SMALL;
  double      initial_dt = 0.0;
  double      max_t      = 0.0;
  bool        save       = false;
  std::string output     = "";
  bool        parallel   = false;
  readArgs(argc, argv, path, name, size, initial_dt, max_t, save, output, parallel);
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Load the model and calculate the trajectory */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  Model* model = new Model(path, name, size, parallel);
  model->compute_local_optimum_for_all_conditions(initial_dt, max_t, save, output);
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Free memory and exit                        */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  delete model;
  model       = NULL;
  clock_t end = clock();
  double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC;
  std::cout << "Elapsed time: " << elapsed_secs << " seconds" << std::endl;
  return EXIT_SUCCESS;
}

/**
 * \brief    Read command line arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    std::string &path
 * \param    std::string &name
 * \param    msize &size
 * \param    double &initial_dt
 * \param    double &max_t
 * \param    bool &save
 * \param    std::string &output
 * \param    bool &parallel
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string &path, std::string &name, msize &size, double &initial_dt, double &max_t, bool &save, std::string &output, bool &parallel )
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
    else if (strcmp(argv[i], "-size") == 0 || strcmp(argv[i], "--model-size") == 0)
    {
      if (i+1 == argc)
      {
        throw std::invalid_argument("> model size value is missing");
      }
      else
      {
        if (strcmp(argv[i+1], "SMALL") == 0 || strcmp(argv[i+1], "small") == 0)
        {
          size = SMALL;
        }
        else if (strcmp(argv[i+1], "GENOME_SCALE") == 0 || strcmp(argv[i+1], "genome_scale") == 0)
        {
          size = GENOME_SCALE;
        }
        else
        {
          throw std::invalid_argument("> wrong value for parameter -size (--model-size).");
        }
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
        else if (strcmp(argv[i], "-save") == 0 || strcmp(argv[i], "--save-trajectory") == 0)
    {
      save = true;
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
    else if (strcmp(argv[i], "-parallel") == 0 || strcmp(argv[i], "--parallel-computing") == 0)
    {
      parallel = true;
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
  std::cout << " Copyright (C) 2024                                                  \n";
  std::cout << " Charles Rocabert                                                    \n";
  std::cout << " Web: https://github.com/charlesrocabert/GBA_Evolution_2/            \n";
  std::cout << "                                                                     \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                     \n";
  std::cout << " This is free software, and you are welcome to redistribute it under \n";
  std::cout << " certain conditions; See the GNU General Public License for details  \n";
  std::cout << "*********************************************************************\n";
  std::cout << "\n";
  std::cout << "Usage: calculate_local_optimum -h or --help\n";
  std::cout << "   or: calculate_local_optimum [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -path, --model-path\n";
  std::cout << "        specify the path of the GBA model to be loaded\n";
  std::cout << "  -name, --model-name\n";
  std::cout << "        specify the name of the GBA model to be loaded\n";
  std::cout << "  -size, --model-size\n";
  std::cout << "        specify the size of the GBA model (SMALL / GENOME_SCALE)\n";
  std::cout << "  -dt, --initial-dt\n";
  std::cout << "        specify the initial gradient timestep\n";
  std::cout << "  -maxt, --max-time\n";
  std::cout << "        specify the maximal time\n";
  std::cout << "  -save, --save-trajectory\n";
  std::cout << "        specify if the trajectory should be saved as output files\n";
  std::cout << "  -output, --output-path\n";
  std::cout << "        specify the path of output files\n";
  std::cout << "  -parallel, --parallel-computing\n";
  std::cout << "        specify if the computation should be parallel\n";
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
  std::cout << " Copyright (C) 2024                                                  \n";
  std::cout << " Charles Rocabert                                                    \n";
  std::cout << " Web: https://github.com/charlesrocabert/GBA_Evolution_2/            \n";
  std::cout << "                                                                     \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                     \n";
  std::cout << " This is free software, and you are welcome to redistribute it under \n";
  std::cout << " certain conditions; See the GNU General Public License for details  \n";
  std::cout << "*********************************************************************\n";
  std::cout << "\n";
}


