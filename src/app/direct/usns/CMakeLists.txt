set(main_target "UnsteadyNavierStokesSolver")

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
    TP petsc metis omp direct config 
    mesh solver boundary common 
    fem export
)

target_include_directories(${main_target} 
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/lib/boundary
    ${PROJECT_SOURCE_DIR}/src/lib/config
    ${PROJECT_SOURCE_DIR}/src/lib/direct
    ${PROJECT_SOURCE_DIR}/src/lib/inverse
    ${PROJECT_SOURCE_DIR}/src/lib/mesh
    ${PROJECT_SOURCE_DIR}/src/lib/solver
    ${PROJECT_SOURCE_DIR}/src/lib/import
    ${PROJECT_SOURCE_DIR}/src/lib/export
    ${PROJECT_SOURCE_DIR}/src/lib/common 
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})