set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(mesh STATIC)
target_sources(mesh
    PRIVATE
        Grid.cpp
        Cell.cpp
        Node.cpp
)

target_include_directories(mesh
    PUBLIC  
        ${PROJECT_SOURCE_DIR}/src/lib/common
        ${PROJECT_SOURCE_DIR}/src/lib/config
        ${PROJECT_SOURCE_DIR}/src/lib/boundary
        ${PROJECT_SOURCE_DIR}/src/lib/mesh
        ${PROJECT_SOURCE_DIR}/src/lib/output
        ${PROJECT_SOURCE_DIR}/src/lib/solver
)

target_compile_options(mesh PUBLIC ${compile_options_tmp})