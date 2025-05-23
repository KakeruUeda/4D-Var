set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(inverse STATIC)
target_sources(inverse
  PRIVATE
    InverseProblem.cpp
    Initialize.cpp
    Adjoint.cpp
    AssemblyAdjoint.cpp
    OptimalCondition.cpp
    Output.cpp
)

target_include_directories(inverse
  PUBLIC
    ${PROJECT_SOURCE_DIR}/src/lib/config
    ${PROJECT_SOURCE_DIR}/src/lib/direct
    ${PROJECT_SOURCE_DIR}/src/lib/mesh
    ${PROJECT_SOURCE_DIR}/src/lib/common
    ${PROJECT_SOURCE_DIR}/src/lib/boundary
    ${PROJECT_SOURCE_DIR}/src/lib/import
    ${PROJECT_SOURCE_DIR}/src/lib/export
    ${PROJECT_SOURCE_DIR}/src/lib/solver
    ${PROJECT_SOURCE_DIR}/src/lib/fem
)

target_compile_options(inverse PUBLIC ${compile_options_tmp})