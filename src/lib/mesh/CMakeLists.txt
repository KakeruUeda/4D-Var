set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(mesh STATIC)
target_sources(mesh
    PRIVATE
        Cell.cpp
)

target_include_directories(mesh 
    PUBLIC  
        ${PROJECT_SOURCE_DIR}/src/lib/tools
        ${PROJECT_SOURCE_DIR}/src/lib/config
)

target_compile_options(mesh PUBLIC ${compile_options_tmp})