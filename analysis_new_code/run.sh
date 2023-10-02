gfortran -fcheck=all -o 1_histogram.x 1_histogram.f
./1_histogram.x -T0 300.d0 -T 1200.d0 -tmin 0 -tmax YYYY -grid -3.14d0 3.14d0 0.1281d0 2.d0 14.d0 0.2448d0 2.d0 14.d0 0.2448d0 -pfrqMD 10 

