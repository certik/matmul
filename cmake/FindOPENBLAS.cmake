find_path(OPENBLAS_INCLUDE_DIR NAMES cblas.h PATHS /usr/include/openblas)

find_library(OPENBLAS_LIBRARY NAMES openblas)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OPENBLAS DEFAULT_MSG OPENBLAS_INCLUDE_DIR
    OPENBLAS_LIBRARY)

add_library(p::openblas INTERFACE IMPORTED)
set_property(TARGET p::openblas PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    ${OPENBLAS_INCLUDE_DIR})
set_property(TARGET p::openblas PROPERTY INTERFACE_LINK_LIBRARIES
    ${OPENBLAS_LIBRARY})
