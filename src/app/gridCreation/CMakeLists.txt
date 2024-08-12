set(main_target "GridCreation")

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
    boundary io common preprocess 
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
    ${PROJECT_SOURCE_DIR}/src/lib/fileIO
    ${PROJECT_SOURCE_DIR}/src/lib/common 
    ${PROJECT_SOURCE_DIR}/src/lib/preprocess
    ${PROJECT_SOURCE_DIR}/src/lib/postprocess 
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})