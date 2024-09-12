mkdir GBAcpp
cp -r ./build ./GBAcpp/.
cp -r ./cmake ./GBAcpp/.
cp -r ./src ./GBAcpp/.
cp -r ./output ./GBAcpp/.
cp -r ./csv_models ./GBAcpp/.
cp -r ./CMakeLists.txt ./GBAcpp/.
zip -r GBAcpp.zip GBAcpp
scp GBAcpp.zip dam82xot@storage.hpc.rz.uni-duesseldorf.de:/gpfs/project/dam82xot/.
rm GBAcpp.zip
rm -rf GBAcpp
