set(main_target "4DVar")

set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_executable(${main_target} main.cpp)

target_sources(${main_target}
    PRIVATE
        main.cpp
)

target_compile_options(${main_target} PUBLIC ${compile_options_tmp})


target_link_libraries( ${main_target} 
    PRIVATE
        TP petsc metis omp direct
        inverse config mesh solver
        boundary output common
)

target_include_directories(${main_target} 
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src/lib/boundary
        ${PROJECT_SOURCE_DIR}/src/lib/config
        ${PROJECT_SOURCE_DIR}/src/lib/direct
        ${PROJECT_SOURCE_DIR}/src/lib/inverse
        ${PROJECT_SOURCE_DIR}/src/lib/mesh
        ${PROJECT_SOURCE_DIR}/src/lib/solver
        ${PROJECT_SOURCE_DIR}/src/lib/output
        ${PROJECT_SOURCE_DIR}/src/lib/common 
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})