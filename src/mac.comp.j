# mac serial
  cd etc
  cp -rp macros.make.macgfortran  macros.make.all  ../
  cd -
  ./makeswm arch=macgfortran  hw=cpu par=serial threading=no >&  za.macgfortran
#  rm -f  macros.make.macgfortran  macros.make.all
