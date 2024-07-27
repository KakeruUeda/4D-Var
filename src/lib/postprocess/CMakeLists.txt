set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(postinverse STATIC)
target_sources(postinverse
    PRIVATE
    PostInverseProblem.cpp
)

target_include_directories(postinverse
    PUBLIC  
        ${PROJECT_SOURCE_DIR}/src/lib/common
        ${PROJECT_SOURCE_DIR}/src/lib/config
        ${PROJECT_SOURCE_DIR}/src/lib/mesh
        ${PROJECT_SOURCE_DIR}/src/lib/boundary
        ${PROJECT_SOURCE_DIR}/src/lib/solver
        ${PROJECT_SOURCE_DIR}/src/lib/fileIO
)

target_compile_options(postinverse PUBLIC ${compile_options_tmp})