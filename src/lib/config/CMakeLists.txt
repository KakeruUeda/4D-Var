set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(config STATIC)
target_sources(config
    PRIVATE
        Config.cpp
)

target_include_directories(config 
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src/lib/tools
)

target_compile_options(config PUBLIC ${compile_options_tmp})