include(FindPackageHandleStandardArgs)

set(IPE_ROOT "" CACHE STRING "Optionally, a root directory to look for the IPE library in.")
mark_as_advanced(IPE_ROOT)

find_path(IPE_INCLUDE_DIR NAMES ipelib.h HINTS ${IPE_ROOT}/include)
find_library(IPE_LIBRARY NAMES ipe HINTS ${IPE_ROOT}/lib /Applications/Ipe.app/Contents/Frameworks)
mark_as_advanced(IPE_INCLUDE_DIR IPE_LIBRARY)

find_package_handle_standard_args(IPE DEFAULT_MSG IPE_LIBRARY IPE_INCLUDE_DIR)

if(IPE_FOUND)
	set(IPE_LIBRARIES ${IPE_LIBRARIES} ${IPE_LIBRARY})
	set(IPE_INCLUDE_DIRS ${IPE_INCLUDE_DIR})

	if(NOT TARGET algutil::ipe)
		util_add_imported_library(algutil::ipe ${IPE_LIBRARY} ${IPE_INCLUDE_DIRS})
	endif()
endif()

