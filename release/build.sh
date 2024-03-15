rm -rf bin
mkdir release
cd release

cmake -D CMAKE_INSTALL_PREFIX=/Users/kakeru/Documents/Projects/DataAssimilation/release/bin \
      ..

make && make install