#!/usr/bin/env python3
# coding: utf-8

#***********************************************************************
# Copyright © 2024 Charles Rocabert
# Web: https://github.com/charlesrocabert/GBA_MMSYN
#
# run_GBAcpp_jobs.py
# ------------------
# Run the model for a given list of conditions.
# (HPC SCRIPT)
#***********************************************************************

import os
import sys

### Build the command line ###
def build_command_line( exec_path, model_path, model_name, condition, dt, maxt, output_path ):
    cmd  = exec_path
    cmd += " -path " + model_path
    cmd += " -name " + model_name
    cmd += " -condition " + condition
    cmd += " -dt " + str(dt)
    cmd += " -maxt " + str(maxt)
    cmd += " -output " + output_path
    return cmd

### Build and run the qsub script ###
def build_and_run_qsub_script( exec_path, model_path, model_name, condition, dt, maxt, output_path ):
    cmd = build_command_line(exec_path, model_path, model_name, condition, dt, maxt, output_path)
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
    EXEC_PATH   = "/gpfs/project/dam82xot/GBAcpp/build/bin/compute_gradient_ascent"
    MODEL_PATH  = "/gpfs/project/dam82xot/GBAcpp/csv_models"
    MODEL_NAMES = ["MMSYN_0000", "MMSYN_0010", "MMSYN_0100", "MMSYN_0110",
                   "MMSYN_1000", "MMSYN_1001", "MMSYN_1010", "MMSYN_1011",
                   "MMSYN_1100", "MMSYN_1101", "MMSYN_1110", "MMSYN_1111"]
    CONDITIONS  = range(1, 31)
    CONDITIONS  = [str(c) for c in CONDITIONS]
    DT          = 0.01
    MAXT        = 1000000.0
    OUTPUT_PATH = "/gpfs/project/dam82xot/GBAcpp/output"
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Run a job for each condition #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    for model_name in MODEL_NAMES:
        for condition in CONDITIONS:
            build_and_run_qsub_script(EXEC_PATH, MODEL_PATH, model_name, condition, DT, MAXT, OUTPUT_PATH)

