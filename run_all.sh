cd cmake
sh make.sh
cd ..
./build/bin/compute_local_optimums -path ./csv_models -name A -dt 0.01 -maxt 10000 -output ./output
./build/bin/compute_local_optimums -path ./csv_models -name B -dt 0.01 -maxt 10000 -output ./output
./build/bin/compute_local_optimums -path ./csv_models -name C -dt 0.01 -maxt 10000 -output ./output
./build/bin/compute_local_optimums -path ./csv_models -name D -dt 0.01 -maxt 10000 -output ./output
./build/bin/compute_local_optimums -path ./csv_models -name EC12b -dt 0.01 -maxt 10000 -output ./output
