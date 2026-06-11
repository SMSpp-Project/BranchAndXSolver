# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ BranchAndXSolver                           #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#                              Donato Meoli                                   #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.21)
cmake_policy(VERSION 3.15)
project(BranchAndXSolver
        VERSION 0.1.0
        DESCRIPTION "SMS++ BranchAndXSolver module"
        HOMEPAGE_URL https://github.com/Fiden99/testBnBKnapsack
        LANGUAGES CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})
set(modNamespace "SMS++")

# Find out if it's being called by the umbrella.
get_directory_property(hasParent PARENT_DIRECTORY)

# This adds the cmake directory to the module search paths,
# allowing us to use our modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Settings ------------------------------------------------------------ #
# Sets the default build type (if none was specified).
# See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
include(BuildType)

# ----- Requirements -------------------------------------------------------- #
# If it's not being called by the umbrella, we need to
# look for the system-installed libraries.
if (NOT hasParent)
    find_package(SMS++ REQUIRED)
    # BinaryKnapsackBlock is NOT required: it only gates the BKB-specific
    # GreedyRelaxationSolver (see below)
    find_package(BinaryKnapsackBlock QUIET)
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_17)

if (MSVC AND BUILD_SHARED_LIBS)
    target_link_options(${modName} INTERFACE "/INCLUDE:SMSpp_force_load_${modName}")
endif ()

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# Note: do not GLOB files here.
target_sources(${modName} PRIVATE
        src/BranchAndXSolver.cpp)

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake).
target_include_directories(
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

# When linking other targets to the library with target_link_libraries(),
# PUBLIC means they will be linked also to the targets that depend on this
# library.
target_link_libraries(${modName} PUBLIC ${modNamespace}::SMS++)

# The BKB-specific GreedyRelaxationSolver is built only when the
# BinaryKnapsackBlock module is part of the build: the BranchAndXSolver
# itself is fully generic and does not require it. (Note that with the
# makefiles the BinaryKnapsackBlock chain is always pulled in, see
# makefile-c / makefile-s.)
if (TARGET ${modNamespace}::BinaryKnapsackBlock)
    message(STATUS "${modName}: BinaryKnapsackBlock found, "
                   "building GreedyRelaxationSolver")
    target_sources(${modName} PRIVATE src/GreedyRelaxationSolver.cpp)
    target_link_libraries(${modName} PUBLIC
            ${modNamespace}::BinaryKnapsackBlock)
    set(BAXSLV_WITH_BKB ON)
else ()
    set(BAXSLV_WITH_BKB OFF)
endif ()

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(${modNamespace}::${modName} ALIAS ${modName})

# ----- Test ---------------------------------------------------------------- #
add_executable(${modName}_test test/test.cpp)
target_link_libraries(${modName}_test PRIVATE ${modNamespace}::${modName})

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS ${modName}
        EXPORT ${modName}Targets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # .dll
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # .so .dylib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}) # .lib

# Install the headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT ${modName}Targets
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${modName}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# --------------------------------------------------------------------------- #
