#!/bin/bash
for j in `seq 1000 500 10000`
#for j in `seq 500 500 500`
do

for i in `cat doc.dat`
 
do
  echo -e "\e[0;31m======================================================================================================================\e[0m"
  echo -e "\e[0;31m                                       ENTERING LAMBDA WINDOW " ${i} Please WAIT...."\e[0m"
  echo -e "\e[0;31m======================================================================================================================\e[0m"
#  mkdir analysis$i
#  cp ../out_file/$i/ti001.out analysis$i
  cd analysis$i
#  rm COLVAR_1
#  cp ../../data$i/COLVAR COLVAR_1
#  sed -i '/^#/d' COLVAR_1
#  sed -n '1p' COLVAR_1 > new && sed -e '1d' COLVAR_1 | sed -n '0~1000p' >>new
#  mv new COLVAR_1
  cp ../analysis_new_code/* . 
#  python extract_dhdl.py >out
#  sed -e '5001d;5002d;10003d;10004d' out >out1
#  mv out1 ti001.out
#  rm out
  cp run_conv.sh run.sh
  sed -i -e "s/YYYY/$j/g" run.sh
  sed -i "s/XXXX/$i/g" 2_dhdl.f90 
  sed -i -e "s/YYYY/$j/g" 2_input.dat
  sh run.sh    # it will generate output file "Pu.dat"
  gfortran 2_dhdl.f90 -o 2_dhdl.x
  ./2_dhdl.x >>../int_dhdl.dat
  cd ../
  echo "            lambda $i done       "
done  
   cp deltaf_conv.f90 deltaf.f90
   sed -i -e "s/YYYY/$j/g" deltaf.f90
   gfortran deltaf.f90
   ./a.out >>out
   mv int_dhdl.dat int_dhdl$j.dat
   mv deltaf.dat deltaf$j.dat
done
