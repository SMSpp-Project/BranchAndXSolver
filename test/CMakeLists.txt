# --------------------------------------------------------------------------- #
#    CMake file for BranchAndXSolver/test                                     #
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
# what its module asks of the linker [see the CMakeLists.txt of the
# module], while an archive has to be taken whole.
if (NOT BUILD_SHARED_LIBS)
    if (MSVC)
        add_link_options("/WHOLEARCHIVE" "/FORCE:MULTIPLE")
    else () # Unix
        if (APPLE)
            add_link_options("-Wl,-all_load")
        else ()
            add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
        endif ()
    endif ()
endif ()

# ----- BranchAndXSolver_unit_test ------------------------------------------ #
add_executable(${modName}_unit_test test.cpp)
target_link_libraries(${modName}_unit_test PRIVATE SMS++::BranchAndXSolver)

if (WIN32 AND BUILD_SHARED_LIBS AND NOT DEFINED ENV{CONDA_BUILD})
    # On Windows with shared libraries enabled copies all runtime DLL
    # dependencies next to the executable in the build directory
    # (POST_BUILD), using CMake's TARGET_RUNTIME_DLLS generator expression
    add_custom_command(TARGET ${modName}_unit_test POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:${modName}_unit_test>
            $<TARGET_FILE_DIR:${modName}_unit_test>
            COMMAND_EXPAND_LISTS
            COMMENT "Copying runtime DLLs next to ${modName}_unit_test.exe")
endif ()

# The tester builds its instances in memory, and reads the BlockSolverConfig
# of the inner Solver out of this directory.
add_test(NAME ${modName}_unit_test
         COMMAND ${modName}_unit_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- the label of this module's tests ------------------------------------ #

# The label lets the module CI select exactly the tests exercising this module
# with `ctest -L <module>` (see tests/cmake/TestLabels.cmake in the umbrella
# for the system-tests counterpart): without it these tests are built and
# never run.
get_property(_tests DIRECTORY PROPERTY TESTS)
if (_tests)
    set_tests_properties(${_tests} PROPERTIES LABELS "${modName}")
endif ()

# --------------------------------------------------------------------------- #
