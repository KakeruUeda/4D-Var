set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(common STATIC)
target_sources(common PRIVATE
    MathTool.cpp
    Spline.cpp
    Polynomial.cpp
    SignalHnadler.cpp
)

target_compile_options(common PUBLIC ${compile_options_tmp})