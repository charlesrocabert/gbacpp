rm -rf ./output/*
scp -r dam82xot@storage.hpc.rz.uni-duesseldorf.de:/gpfs/project/dam82xot/GBAcpp/output/* ./output/.

Rscript plot_trajectory.R 1
Rscript plot_trajectory.R 2
Rscript plot_trajectory.R 3
Rscript plot_trajectory.R 4
Rscript plot_trajectory.R 5

