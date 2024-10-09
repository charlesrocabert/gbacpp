#!/usr/bin/env python3
# coding: utf-8

MODELS      = ["MMSYN_0000", "MMSYN_1000", "MMSYN_0100", "MMSYN_0010",
               "MMSYN_0001", "MMSYN_1100", "MMSYN_1010", "MMSYN_1001",
               "MMSYN_0110", "MMSYN_0101", "MMSYN_0011", "MMSYN_1110",
               "MMSYN_1101", "MMSYN_1011", "MMSYN_0111", "MMSYN_1111"]
CONDITIONS  = range(1, 31)
CONDITIONS  = [str(c) for c in CONDITIONS]
LABELS      = ["b", "c", "f", "p", "state", "v"]

outputs = {}
for model in MODELS:
    for label in LABELS:
        filename                 = "./output/"+model+"_"+label+"_optimum.csv"
        outputs[model+"_"+label] = open(filename, "w")

for model in MODELS:
    for label in LABELS:
        header_written = False
        for condition in CONDITIONS:
            filename = "./hpc_download/"+model+"_"+condition+"_"+label+"_optimum.csv"
            try:
                f = open(filename, "r")
                if not header_written:
                    outputs[model+"_"+label].write(f.readline())
                    header_written = True
                else:
                    f.readline()
                outputs[model+"_"+label].write(f.readline())
                f.close()
            except:
                print("> Warning: file "+filename+" not found")

for item in outputs.items():
    item[1].close()