if (SLAP_CLANG_TIDY AND (CMAKE_C_COMPILER_ID STREQUAL "Clang"))
  message(STATUS "Running clang-tidy on source directory.")
  set(CMAKE_C_CLANG_TIDY
    clang-tidy;
    -export-fixes=tidy-fixes.yaml;
    -header-filter=.;
    --config-file=${PROJECT_SOURCE_DIR}/.clang-tidy;
    )
endif ()


# include(FindCPM)
# CPMAddPackage(
#   NAME Eigen
#   VERSION 3.4.0
#   URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
#   # Eigen's CMakelists are not intended for library use
#   DOWNLOAD_ONLY YES
# )
find_package(Eigen3)
# if(Eigen_ADDED)
#   add_library(Eigen::Eigen INTERFACE IMPORTED)
#   target_include_directories(Eigen::Eigen INTERFACE ${Eigen_SOURCE_DIR})
# endif()

add_library(slap
  matrix.h
  matrix.c

  errors.h
  errors.c

  matrix_checks.h
  matrix_checks.c

  printing.h
  printing.c

  new_matrix.h
  new_matrix.c

  copy_matrix.h
  copy_matrix.c

  unary_ops.h
  unary_ops.c

  binary_ops.h
  binary_ops.c

  vector_ops.h
  vector_ops.c

  vector_products.h
  vector_products.c

  iterator.h
  iterator.c

  strided_matrix.h
  strided_matrix.c

  function_mapping.h
  function_mapping.c

  matmul.h
  matmul.c
  matmul_eigen.cpp

  cholesky.h
  cholesky.c qr.c qr.h tri.c tri.h)

target_include_directories(slap
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:"include">
)

target_compile_definitions(slap PUBLIC SLAP_FLOAT=${SLAP_FLOAT})

# Link math library
if (NOT APPLE AND NOT WIN32)
  target_link_libraries(slap PUBLIC m)
  target_link_libraries(slap PRIVATE Eigen3::Eigen)
endif()

add_library(slap::slap ALIAS slap)