set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(direct STATIC)
target_sources(direct
    PRIVATE
        DirectProblem.cpp
        MatAssySNS.cpp
        MatAssyUSNS.cpp
)

target_include_directories(direct
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src/lib/config
        ${PROJECT_SOURCE_DIR}/src/lib/mesh
        ${PROJECT_SOURCE_DIR}/src/lib/common
)

target_compile_options(direct PUBLIC ${compile_options_tmp})