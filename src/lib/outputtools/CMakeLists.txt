set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(outputtools STATIC)
target_sources(outputtools
    PRIVATE
        VTIStructuredGrid.cpp
)

target_include_directories(outputtools
    PUBLIC  
        ${PROJECT_SOURCE_DIR}/src/lib/common
        ${PROJECT_SOURCE_DIR}/src/lib/config
        ${PROJECT_SOURCE_DIR}/src/lib/boundary
)

target_compile_options(outputtools PUBLIC ${compile_options_tmp})