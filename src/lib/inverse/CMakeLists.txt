set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(inverse STATIC)
target_sources(inverse
    PRIVATE
        InverseProblem.cpp
        Initialize.cpp
        Adjoint.cpp
        AssemblyAdjoint.cpp
)

target_include_directories(inverse
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src/lib/config
        ${PROJECT_SOURCE_DIR}/src/lib/direct
        ${PROJECT_SOURCE_DIR}/src/lib/mesh
        ${PROJECT_SOURCE_DIR}/src/lib/common
        ${PROJECT_SOURCE_DIR}/src/lib/boundary
        ${PROJECT_SOURCE_DIR}/src/lib/fileIO
        ${PROJECT_SOURCE_DIR}/src/lib/solver
)

target_compile_options(inverse PUBLIC ${compile_options_tmp})