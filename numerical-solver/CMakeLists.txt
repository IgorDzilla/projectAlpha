cmake_minimum_required(VERSION 3.16)
project(numerical_solver LANGUAGES CXX)

add_library(numerics STATIC
  src/problems/RiemannProblem1D.cpp
  src/methods/FiniteVolume.cpp
  src/io/CsvResultWriter.cpp
  src/models/EulerEquations.cpp
  src/methods/riemann/ExactEulerRiemannSolver.cpp
)
target_include_directories(numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(numerics PUBLIC cxx_std_17)

set(CONFIG_READER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../utils")
add_library(config-reader STATIC
  "${CONFIG_READER_DIR}/toml-processor.cpp"
  "${CONFIG_READER_DIR}/math-expression.cpp"
)
target_include_directories(config-reader PUBLIC "${CONFIG_READER_DIR}")
target_compile_features(config-reader PUBLIC cxx_std_17)

# TOML остаётся в слое приложения; старый wave_1d reader не меняется.
add_library(euler-app STATIC
  apps/solver-cli/config/EulerConfig.cpp
  apps/solver-cli/EulerRun.cpp
)
target_include_directories(euler-app PUBLIC apps/solver-cli "${CONFIG_READER_DIR}")
target_link_libraries(euler-app PUBLIC numerics)
add_executable(solver-cli apps/solver-cli/main.cpp)
target_link_libraries(solver-cli PRIVATE euler-app)

include(CTest)
if(BUILD_TESTING)
  add_executable(finite-volume-test tests/finite-volume-test.cpp)
  target_link_libraries(finite-volume-test PRIVATE numerics)
  add_test(NAME finite-volume COMMAND finite-volume-test)

  add_executable(euler-config-test tests/euler-config-test.cpp)
  target_link_libraries(euler-config-test PRIVATE euler-app)
  add_test(NAME euler-config COMMAND euler-config-test
    "${CMAKE_CURRENT_SOURCE_DIR}/../configs" "${CMAKE_CURRENT_BINARY_DIR}/config-test-results")

  add_executable(check-euler-output tests/check-euler-output.cpp)
  target_link_libraries(check-euler-output PRIVATE euler-app)
  foreach(case RANGE 1 5)
    add_test(NAME cli-toro-${case} COMMAND ${CMAKE_COMMAND}
      "-DCLI=$<TARGET_FILE:solver-cli>" "-DCHECKER=$<TARGET_FILE:check-euler-output>"
      "-DCONFIG=${CMAKE_CURRENT_SOURCE_DIR}/../configs/toro_${case}.toml"
      "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/cli-test-results/toro_${case}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cli-test.cmake")
  endforeach()

  add_executable(core-test tests/core-test.cpp)
  target_link_libraries(core-test PRIVATE numerics)
  add_test(NAME core COMMAND core-test)

  add_executable(solver-test tests/solver-test.cpp)
  target_link_libraries(solver-test PRIVATE numerics)
  add_test(NAME solver COMMAND solver-test)

  add_executable(riemann-toro-test tests/riemann-toro-test.cpp)
  target_link_libraries(riemann-toro-test PRIVATE numerics)
  add_test(NAME riemann-toro COMMAND riemann-toro-test)

  add_executable(riemann-test tests/riemann-test.cpp)
  target_link_libraries(riemann-test PRIVATE numerics)
  add_test(NAME riemann COMMAND riemann-test)

  add_executable(config-reader-test "${CONFIG_READER_DIR}/tests/toml-processor-test.cpp")
  target_link_libraries(config-reader-test PRIVATE config-reader)
  add_test(NAME config-reader COMMAND config-reader-test)
  set_tests_properties(config-reader PROPERTIES WORKING_DIRECTORY "${CONFIG_READER_DIR}")
endif()
