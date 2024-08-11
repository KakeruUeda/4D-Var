set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(fem STATIC)
target_sources(fem
  PRIVATE
  FEM.cpp
)

target_include_directories(fem
  PUBLIC  
    ${PROJECT_SOURCE_DIR}/src/lib/common
    ${PROJECT_SOURCE_DIR}/src/lib/config
    ${PROJECT_SOURCE_DIR}/src/lib/boundary
    ${PROJECT_SOURCE_DIR}/src/lib/mesh
    ${PROJECT_SOURCE_DIR}/src/lib/fileIO
    ${PROJECT_SOURCE_DIR}/src/lib/import
    ${PROJECT_SOURCE_DIR}/src/lib/solver
    ${PROJECT_SOURCE_DIR}/src/lib/preprocess
    ${PROJECT_SOURCE_DIR}/src/lib/postprocess
)

target_compile_options(fem PUBLIC ${compile_options_tmp})