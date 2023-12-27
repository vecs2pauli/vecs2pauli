cmake_minimum_required(VERSION 3.9)

# Set a name and a version number for your project:
project(vecs2pauli VERSION 0.0.1 LANGUAGES CXX)

# Initialize some default paths
include(GNUInstallDirs)


# Define the minimum C++ standard that is required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Compilation options
option(BUILD_PYTHON "Enable building of Python bindings" ON)

# compile the library
add_subdirectory(src)


# compile the application
add_subdirectory(app)


# include binarylinalg library
add_subdirectory(extern/binarylinalg)



# compile the tests
#include(CTest)
#if(BUILD_TESTING)
#  find_package(Catch2 REQUIRED)
#  include(Catch)
#  add_subdirectory(tests)
#endif()


if(BUILD_PYTHON)
  # Add Python bindings
  find_package(pybind11)
  # Compile the Pybind11 module
  pybind11_add_module(_vecs2pauli python/vecs2pauli/_vecs2pauli.cpp)
  target_link_libraries(_vecs2pauli PRIVATE vecs2pauli binarylinalg)
  target_include_directories(_vecs2pauli PRIVATE "${PROJECT_SOURCE_DIR}/extern/")

  # Install the Python module shared library
  install(TARGETS _vecs2pauli DESTINATION .)
endif()

# Add an alias target for use if this project is included as a subproject in another project
add_library(vecs2pauli::vecs2pauli ALIAS vecs2pauli)

# Install targets and configuration
install(
  TARGETS vecs2pauli
  EXPORT vecs2pauli-config
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(
  EXPORT vecs2pauli-config
  NAMESPACE vecs2pauli::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/vecs2pauli
)

install(
  DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# This prints a summary of found dependencies
include(FeatureSummary)
feature_summary(WHAT ALL)



