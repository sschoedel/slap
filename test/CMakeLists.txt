# function add_ulqr_test(name)
#
# Adds a new executable test called <name>_test.
# Assumes thes source code is in a file called <name>_test.
function (add_slap_test name)
  set(TEST_NAME ${name}_test)
  add_executable(${TEST_NAME}
    ${TEST_NAME}.cpp
  )
  target_link_libraries(${TEST_NAME}
    PRIVATE 
    slap
    gtest::gtest
    Eigen3::Eigen
  )
  if (NOT APPLE AND NOT WIN32)
    target_link_libraries(${TEST_NAME} PUBLIC m)
  endif()
  gtest_discover_tests(${TEST_NAME})
endfunction()

find_package(Eigen3)

# Add tests
add_slap_test(eigen_backend)
add_slap_test(matrix)
add_slap_test(vector)
add_slap_test(submatrix)
add_slap_test(linear_algebra)
add_slap_test(errors)