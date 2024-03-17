set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(postprocess STATIC)
target_sources(postprocess
    PRIVATE
        PostDirect.cpp
        PostInverse.cpp
)

target_include_directories(postprocess 
    PUBLIC
)

target_compile_options(postprocess PUBLIC ${compile_options_tmp})