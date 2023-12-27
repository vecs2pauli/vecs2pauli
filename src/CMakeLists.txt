add_library(vecs2pauli vecs2pauli.cpp)
target_include_directories(vecs2pauli PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
