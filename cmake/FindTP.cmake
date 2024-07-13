# https://qiita.com/shohirose/items/d9bda00a39a113965c5c

find_path(TP_INCLUDE_DIR TextParser.h
  PATHS
    ${TP_DIR}
    $ENV{TP_DIR}
    /usr
    /usr/local
  PATH_SUFFIXES
    include
)

find_library(TP_LIBRARY 
NAMES TP
PATHS
    ${TP_DIR}
    $ENV{TP_DIR}
    /usr
    /usr/local
  PATH_SUFFIXES
    lib
) 

mark_as_advanced(
  TP_INCLUDE_DIR
  TP_LIBRARY     # ヘッダーのみのライブラリの場合は不要
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TP
  REQUIRED_VARS
    TP_INCLUDE_DIR
    TP_LIBRARY      # ヘッダーのみのライブラリの場合は不要
  )

if(TP_FOUND AND NOT TARGET TP)
  add_library(TP UNKNOWN IMPORTED)
  set_target_properties(TP PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES ["C"|"CXX"]  # ヘッダーのみのライブラリの場合は不要
    IMPORTED_LOCATION "${TP_LIBRARY}"       # ヘッダーのみのライブラリの場合は不要
    INTERFACE_INCLUDE_DIRECTORIES "${TP_INCLUDE_DIR}"
    )
endif()