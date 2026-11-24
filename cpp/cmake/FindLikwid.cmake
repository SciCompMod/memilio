# Finds the likwid library for performance measurements with markers, see https://github.com/RRZE-HPC/likwid
#
# Defines the imported target likwid::likwid. Set Likwid_ROOT, LIKWID_ROOT or LIKWID_DIR to the likwid installation if it
# is not found automatically. Directories in LIBRARY_PATH and CPATH (e.g. set by environment modules) are searched as
# well.

find_library(LIKWID_LIBRARY likwid
    HINTS ENV LIKWID_ROOT ENV LIKWID_DIR ENV LIBRARY_PATH
    PATH_SUFFIXES lib lib64
)
find_path(LIKWID_INCLUDE_DIR likwid-marker.h
    HINTS ENV LIKWID_ROOT ENV LIKWID_DIR ENV CPATH
    PATH_SUFFIXES include
)
mark_as_advanced(LIKWID_LIBRARY LIKWID_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Likwid REQUIRED_VARS LIKWID_LIBRARY LIKWID_INCLUDE_DIR)

if(Likwid_FOUND AND NOT TARGET likwid::likwid)
    add_library(likwid::likwid UNKNOWN IMPORTED)
    set_target_properties(likwid::likwid PROPERTIES
        IMPORTED_LOCATION "${LIKWID_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${LIKWID_INCLUDE_DIR}"
    )
endif()
