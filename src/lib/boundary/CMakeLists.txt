set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(boundary STATIC)
target_sources(boundary
    PRIVATE
    Boundary.cpp
)

target_include_directories(boundary
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src/lib/config
        ${PROJECT_SOURCE_DIR}/src/lib/mesh
        ${PROJECT_SOURCE_DIR}/src/lib/common
        ${PROJECT_SOURCE_DIR}/src/lib/solver
)

target_compile_options(boundary PUBLIC ${compile_options_tmp})