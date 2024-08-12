set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(export STATIC)
target_sources(export
  PRIVATE
    ExportDAT.cpp
    ExportVTI.cpp
    ExportVTU.cpp
)

target_include_directories(export
  PUBLIC  
    ${PROJECT_SOURCE_DIR}/src/lib/common
    ${PROJECT_SOURCE_DIR}/src/lib/config
    ${PROJECT_SOURCE_DIR}/src/lib/boundary
    ${PROJECT_SOURCE_DIR}/src/lib/mesh
    ${PROJECT_SOURCE_DIR}/src/lib/solver
    ${PROJECT_SOURCE_DIR}/src/lib/fem
    ${PROJECT_SOURCE_DIR}/src/lib/import
)

target_compile_options(export PUBLIC ${compile_options_tmp})