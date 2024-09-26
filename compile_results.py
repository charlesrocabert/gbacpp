#!/usr/bin/env python3
# coding: utf-8

MODEL       = "newextended16"
CONDITIONS  = range(1, 201)
CONDITIONS  = [str(c) for c in CONDITIONS]
LABELS      = ["b", "c", "f", "p", "state", "v"]

outputs = {}
for label in LABELS:
    filename       = "./output/"+MODEL+"_"+label+"_optimum.csv"
    outputs[label] = open(filename, "w")

for condition in CONDITIONS:
    for label in LABELS:
        filename = "./hpc_download/"+MODEL+"_"+condition+"_"+label+"_optimum.csv"
        f = open(filename, "r")
        if condition == "1":
            outputs[label].write(f.readline())
        else:
            f.readline()
        outputs[label].write(f.readline())
        f.close()

for label in LABELS:
    outputs[label].close()