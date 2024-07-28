set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(io STATIC)
target_sources(io
    PRIVATE
        VTI.cpp
        VTU.cpp
        BIN.cpp
)

target_include_directories(io
    PUBLIC  
        ${PROJECT_SOURCE_DIR}/src/lib/common
        ${PROJECT_SOURCE_DIR}/src/lib/config
        ${PROJECT_SOURCE_DIR}/src/lib/boundary
        ${PROJECT_SOURCE_DIR}/src/lib/mesh
        ${PROJECT_SOURCE_DIR}/src/lib/solver
)

target_compile_options(io PUBLIC ${compile_options_tmp})
