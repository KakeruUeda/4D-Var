set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(mesh STATIC)
target_sources(mesh
    PRIVATE
        Cell.cpp
        Grid.cpp
)

target_include_directories(mesh
    PUBLIC  
        ${PROJECT_SOURCE_DIR}/src/lib/common
        ${PROJECT_SOURCE_DIR}/src/lib/config
        ${PROJECT_SOURCE_DIR}/src/lib/boundary
)

target_compile_options(mesh PUBLIC ${compile_options_tmp})