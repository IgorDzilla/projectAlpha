# OUTPUT — выделенный каталог внутри дерева сборки, только для этого теста.
file(REMOVE_RECURSE "${OUTPUT}")
execute_process(COMMAND "${CLI}" "${CONFIG}" --output-dir "${OUTPUT}"
  RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "CLI failed: ${stdout}\n${stderr}")
endif()
message(STATUS "${stdout}")
execute_process(COMMAND "${CHECKER}" "${CONFIG}" "${OUTPUT}"
  RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "CSV check failed: ${stdout}\n${stderr}")
endif()
message(STATUS "${stdout}")
