rm -rf bin
mkdir build
cd build

cmake -D CMAKE_INSTALL_PREFIX=/Users/kakeru/Documents/Projects/DataAssimilation/bin \
      ..

make && make install