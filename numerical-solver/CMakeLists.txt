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

# Единственная внешняя зависимость газового ввода — заголовок toml.hpp.
set(TOML_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../utils")

# Пользовательская библиотека сохраняется отдельно от газового приложения.
add_library(config-reader STATIC
  "${TOML_INCLUDE_DIR}/toml-processor.cpp"
  "${TOML_INCLUDE_DIR}/math-expression.cpp"
)
target_include_directories(config-reader PUBLIC "${TOML_INCLUDE_DIR}")
target_compile_features(config-reader PUBLIC cxx_std_17)

add_library(euler-app STATIC
  apps/solver-cli/config/EulerConfig.cpp
  apps/solver-cli/EulerRun.cpp
)
target_include_directories(euler-app PUBLIC apps/solver-cli "${TOML_INCLUDE_DIR}")
target_link_libraries(euler-app PUBLIC numerics)
add_executable(solver-cli apps/solver-cli/main.cpp)
target_link_libraries(solver-cli PRIVATE euler-app)

include(CTest)
if(BUILD_TESTING)
  add_executable(finite-volume-test tests/finite-volume-test.cpp)
  target_link_libraries(finite-volume-test PRIVATE euler-app)
  add_test(NAME finite-volume COMMAND finite-volume-test "${CMAKE_CURRENT_SOURCE_DIR}/../configs/toro_1.toml")

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

  add_executable(riemann-test tests/riemann-test.cpp)
  target_link_libraries(riemann-test PRIVATE euler-app)
  add_test(NAME riemann COMMAND riemann-test
    "${CMAKE_CURRENT_SOURCE_DIR}/../configs/toro_1.toml"
    "${CMAKE_CURRENT_SOURCE_DIR}/../configs/toro_3.toml")

  add_executable(config-reader-test "${TOML_INCLUDE_DIR}/tests/toml-processor-test.cpp")
  target_link_libraries(config-reader-test PRIVATE config-reader)
  add_test(NAME config-reader COMMAND config-reader-test)
  set_tests_properties(config-reader PROPERTIES WORKING_DIRECTORY "${TOML_INCLUDE_DIR}")
endif()
