export CXX=$(ls $CONDA_PREFIX/bin/*g++* $CONDA_PREFIX/bin/*clang++* 2>/dev/null | head -1)
export CXX14=$CXX
export CXX17=$CXX
# add to activate.sh
mkdir -p $CONDA_PREFIX/lib/R/etc/
echo "CXX = $CXX" > $CONDA_PREFIX/lib/R/etc/Makevars.site
echo "CXX14 = $CXX" >> $CONDA_PREFIX/lib/R/etc/Makevars.site
echo "CXX17 = $CXX" >> $CONDA_PREFIX/lib/R/etc/Makevars.site
