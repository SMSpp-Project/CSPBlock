# --------------------------------------------------------------------------- #
#    CMake file for CSPBlock/test                                             #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# The classes of the modules register themselves in the factories from
# static objects, which nothing here names: a shared library is kept by
# what its module asks of the linker, while an archive has to be taken
# whole [see SMSppModule.cmake of SMS++].
smspp_keep_modules()

# ----- CSPBlock_unit_test -------------------------------------------------- #

add_executable(${modName}_unit_test test.cpp)
target_link_libraries(${modName}_unit_test PRIVATE ${modNamespace}::${modName})

smspp_copy_runtime(${modName}_unit_test)

add_test(NAME ${modName}_unit_test
         COMMAND ${modName}_unit_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# The label lets the module CI select exactly the tests exercising this
# module with `ctest -L <module>` (see tests/cmake/TestLabels.cmake in the
# umbrella for the system-tests counterpart).
set_tests_properties(${modName}_unit_test PROPERTIES LABELS "${modName}")

# ----- CSPBlock_battery_test ----------------------------------------------- #

# one test per line of the batch file: the family, the size and the number
# of the solutions of the instance [see battery.cpp]
add_executable(${modName}_battery_test battery.cpp)
target_link_libraries(${modName}_battery_test PRIVATE
                      ${modNamespace}::${modName})

smspp_copy_runtime(${modName}_battery_test)

file(STRINGS batch batchLines REGEX "^[^#]")
foreach (line IN LISTS batchLines)
    separate_arguments(args UNIX_COMMAND "${line}")
    list(GET args 0 family)
    list(GET args 1 size)
    add_test(NAME ${modName}_battery_${family}_${size}
             COMMAND ${modName}_battery_test ${args}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${modName}_battery_${family}_${size} PROPERTIES
                         LABELS "${modName}")
endforeach ()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS batch)

# ----- FlatZincCSPSolver_test ---------------------------------------------- #

# the Solver that solves a CSPBlock through an external FlatZinc solver; the
# test is skipped (exit 77) when the backend of the configuration is not
# installed, so it builds and runs without a FlatZinc solver
add_executable(${modName}_FlatZincCSPSolver_test FlatZincCSPSolver_test.cpp)
target_link_libraries(${modName}_FlatZincCSPSolver_test PRIVATE
                      ${modNamespace}::${modName})

smspp_copy_runtime(${modName}_FlatZincCSPSolver_test)

# the configuration files go next to the test, which reads them
foreach (cfg FlatZincCSPSolverBSCfg.txt FlatZincCSPSolverCfg.txt)
    configure_file(${cfg} ${CMAKE_CURRENT_BINARY_DIR}/${cfg} COPYONLY)
endforeach ()

add_test(NAME ${modName}_FlatZincCSPSolver_test
         COMMAND ${modName}_FlatZincCSPSolver_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(${modName}_FlatZincCSPSolver_test PROPERTIES
                     LABELS "${modName}" SKIP_RETURN_CODE 77)

# --------------------------------------------------------------------------- #
