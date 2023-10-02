gfortran -fcheck=all -o 1_histogram.x 1_histogram.f
./1_histogram.x -T0 300.d0 -T 1200.d0 -tmin 0 -tmax YYYY -grid 2.d0 40.d0 0.5d0 2.d0 40.d0 0.5d0 2.d0 40.d0 0.5d0 2.d0 40.d0 0.5d0 -pfrqMD 10 

