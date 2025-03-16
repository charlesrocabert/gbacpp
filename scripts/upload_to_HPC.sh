mkdir GBAcpp
mkdir ./GBAcpp/output
cp -r ../build ./GBAcpp/.
cp -r ../cmake ./GBAcpp/.
cp -r ../src ./GBAcpp/.
cp -r ../csv_models ./GBAcpp/.
cp -r ../CMakeLists.txt ./GBAcpp/.
cp -r run_GBAcpp_jobs.py ./GBAcpp/.
zip -r GBAcpp.zip GBAcpp
scp GBAcpp.zip dam82xot@storage.hpc.rz.uni-duesseldorf.de:/gpfs/project/dam82xot/.
scp run_GBAcpp_jobs.py dam82xot@storage.hpc.rz.uni-duesseldorf.de:./
rm GBAcpp.zip
rm -rf GBAcpp
