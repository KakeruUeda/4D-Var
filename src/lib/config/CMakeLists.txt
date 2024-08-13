set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(config STATIC)
target_sources(config
  PRIVATE
    Config.cpp
    ReadTP.cpp
    StrBoundarySetting.cpp
    StrGridSetting.cpp
    GridFiltering.cpp
    TextReaderInterface.cpp
    TextReaderGridCreation.cpp
    TextReaderUSNS.cpp
)

target_include_directories(config 
  PUBLIC
    ${PROJECT_SOURCE_DIR}/src/lib/import
    ${PROJECT_SOURCE_DIR}/src/lib/common
)

target_compile_options(config PUBLIC ${compile_options_tmp})