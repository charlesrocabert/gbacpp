#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#***********************************************************************
# GBAcpp (Growth Balance Analysis for C++)
# Copyright © 2024-2025 Charles Rocabert
# Web: https://github.com/charlesrocabert/GBAcpp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

"""
Filename: run_GBAcpp_jobs.py
Author: Charles Rocabert
Date: 2024-16-12
Description:
    HPC jobs helper script for the module GBAcpp.
License: GNU General Public License v3.0
Copyright: © 2024-2025 Charles Rocabert
"""

import os
import sys

### Read the list of conditions from CSV file ###
def load_conditions_from_csv( model_path, model_name ):
    f          = open(model_path+"/"+model_name+"/conditions.csv", "r")
    l          = f.readline()
    conditions = l.strip("\n").split(";") 
    conditions.pop(0)
    f.close()
    return conditions

### Build the command line ###
def build_command_line( exec_path, model_path, model_name, condition, dt, maxt, max_mu_count, output_path ):
    cmd  = exec_path
    cmd += " -path " + model_path
    cmd += " -name " + model_name
    cmd += " -condition " + condition
    cmd += " -dt " + str(dt)
    cmd += " -maxt " + str(maxt)
    cmd += " -max-mu-count "+str(max_mu_count)
    #cmd += " -save"
    cmd += " -output " + output_path
    return cmd

### Build and run the qsub script ###
def build_and_run_qsub_script( exec_path, model_path, model_name, condition, dt, maxt, max_mu_count, output_path ):
    cmd = build_command_line(exec_path, model_path, model_name, condition, dt, maxt, max_mu_count, output_path)
    f   = open("GBAcpp_job.sh", "w")
    f.write("#!/bin/bash\n")
    f.write("#PBS -l select=1:ncpus=1:mem=4GB\n")
    f.write("#PBS -l walltime=24:00:00\n")
    f.write("#PBS -A GBA_Evolution\n")
    f.write("\n")
    f.write(cmd)
    f.write("\n")
    f.close()
    print("> Submitting job for condition " + condition)
    print("  --> " + cmd)
    os.system("qsub GBAcpp_job.sh")


##################
#      MAIN      #
##################

if __name__ == "__main__":

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Define main parameters       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    EXEC_PATH    = "/gpfs/project/dam82xot/GBAcpp/build/bin/compute_gradient_ascent"
    MODEL_PATH   = "/gpfs/project/dam82xot/GBAcpp/csv_models"
    MODEL_NAMES  = ["mmsyn_fcr_v2"]
    DT           = 0.01
    MAXT         = 1000000.0
    MAX_MU_COUNT = 1000
    OUTPUT_PATH  = "/gpfs/project/dam82xot/GBAcpp/output"
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Run a job for each condition #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    for model_name in MODEL_NAMES:
        conditions = load_conditions_from_csv(MODEL_PATH, model_name)
        for condition in conditions:
            build_and_run_qsub_script(EXEC_PATH, MODEL_PATH, model_name, condition, DT, MAXT, MAX_MU_COUNT, OUTPUT_PATH)

