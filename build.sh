rm -rf bin
mkdir build
cd build

cmake -D CMAKE_INSTALL_PREFIX=/Users/kakeru/Documents/Projects/lab/4D-Var/bin \
      -D TP_DIR=/usr/local/TextParser \
      -D EIGEN_DIR=/opt/homebrew/include/eigen3 \
      ..

make && make install