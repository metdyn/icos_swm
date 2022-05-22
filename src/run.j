#!/bin/bash
# usage:   
#   compile               ---  sh run.j     
#   grid generation only  ---  sh run.j 1  DIR
#   SW test only          ---  sh run.j 2  DIR
#   grid + SW-test        ---           3  DIR
machine="theia"; post='theiapgi'
machine="sjet" ; post='sjetintel'
machine="theia"; post='theiaintel'
machine="hera";  post='heraintel'
machine="mac"  ; post='macgfortran'



cat > ${machine}.comp.j <<EOF
# ${machine} serial
  cd etc
  cp -rp macros.make.${post}  macros.make.all  ../
  cd -
  ./makeswm arch=$post  hw=cpu par=serial threading=no >&  za.$post
#  rm -f  macros.make.$post  macros.make.all
EOF
if [ $# -eq 0 ]; then
  echo 'compile swm.x'
  sh ${machine}.comp.j
  cat za.*
  exit 0
fi


if [ $machine == "mac" ]; then
  export OMP_NUM_THREADS=1
fi


cur_dir=`pwd`
store_dir="${cur_dir}/../run"
exe_dir="${cur_dir}/../run_${post}_cpu_serial_ompno"
if [ $# -eq 2 ]; then
  test_dir=${exe_dir}/$2
else
  test_dir=${exe_dir}/.
fi
mkdir -p ${test_dir}
#
if [ -e ${test_dir}/SWMnamelist ]; then
  echo "-e exist: ${test_dir}/SWMnamelist"
else
  echo "-e does not exist: ${test_dir}/SWMnamelist"
  echo "copy SWMnamelist from run/."
  cp -rp ${store_dir}/SWMnamelist  ${test_dir}/.
fi



if [ $# -eq 1 ]; then
  echo '\$# -eq 1  is not de-activiated'
  exit
elif [ $# -eq 2 ]; then
  if [ $1 -le 0 ]; then
      echo 'not implemented'
      exit
  elif [ $1 -eq 1 ]; then      # case 1,2,3
      echo ''
      echo "generate grid"
      echo ''
      echo "cd $test_dir"
      cd $test_dir
      mv -f  z.check*0 ~/trash/.	  
      ${exe_dir}/grid   > z.check.grid
  elif [ $1 -eq 2 ]; then
      echo ''
      echo 'SW-test'
      echo "cd $test_dir"
      echo ''
      cd $test_dir
      # mv -f  d.* ~/trash/.
      echo 'execute swa.x'
      ${exe_dir}/swa.x  > out.sw
  elif [ $1 -eq 3 ]; then      # case 1,2,3
      echo ''
      echo "generate grid"
      echo 'execute ./grid'
      echo "cd $test_dir"
      echo ''
      cd $test_dir
      mv -f  z.check*0 ~/trash/.
      ${exe_dir}/grid  > z.check.grid
      #
      echo ''
      echo 'SW-test'
      echo 'execute swa.x'
      echo ''
      # mv -f  d.* ~/trash/.
      ${exe_dir}/swa.x > out.sw
  else
      echo 'wrong input option, STOP'
      exit
  fi
fi
