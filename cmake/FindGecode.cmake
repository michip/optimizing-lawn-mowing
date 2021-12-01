# find a Gecode installation
include(FindPackageHandleStandardArgs)

if(NOT TARGET algutil::Gecode)
find_path(GECODE_INCLUDE_DIR
	NAMES
		gecode/search.hh
	HINTS
		${GECODE_ROOT}/include
		/usr/local/include
		/opt/local/include
)

set(_GECODE_LIBS 
	support
	kernel
	int
	float
	set
	search
	minimodel
	gist
	driver
	flatzinc
)

set(_LIB_VARS "")
set(_LAST_FOUND_LIB_FILE "")
foreach(_LIB IN LISTS _GECODE_LIBS)
	if(NOT "${_LIB}" STREQUAL "")
		set(_LIB_VARS "${_LIB_VARS};_GECODE_${_LIB}_LIBRARY")
		find_library(_GECODE_${_LIB}_LIBRARY
			NAMES
				"gecode${_LIB}"
			HINTS
				${GECODE_ROOT}
				${GECODE_ROOT}/lib
				/opt/local/lib
				/usr/local/lib
		)

		if(_GECODE_${_LIB}_LIBRARY)
			set(_LAST_FOUND_LIB_FILE "${_GECODE_${_LIB}_LIBRARY}")
		endif()
	endif()
endforeach()

if(GECODE_INCLUDE_DIR AND _GECODE_support_LIBRARY AND _GECODE_kernel_LIBRARY)
	set(GECODE_FOUND On)
else()
	set(GECODE_FOUND Off)
endif()

find_package_handle_standard_args(GECODE DEFAULT_MSG GECODE_INCLUDE_DIR _GECODE_support_LIBRARY _GECODE_kernel_LIBRARY)

if(GECODE_FOUND)
	util_add_imported_library(algutil::Gecode::support ${_GECODE_support_LIBRARY} ${GECODE_INCLUDE_DIR})
	util_add_imported_library(algutil::Gecode ${_LAST_FOUND_LIB_FILE} ${GECODE_INCLUDE_DIR})
	util_add_imported_library(algutil::Gecode::kernel ${_GECODE_kernel_LIBRARY} ${GECODE_INCLUDE_DIR})
	util_imported_link_libraries(algutil::Gecode::kernel algutil::Gecode::support)
	util_imported_link_libraries(algutil::Gecode algutil::Gecode::kernel)

	if(_GECODE_int_LIBRARY)
		util_add_imported_library(algutil::Gecode::int ${_GECODE_int_LIBRARY} ${GECODE_INCLUDE_DIR})
		util_imported_link_libraries(algutil::Gecode::int algutil::Gecode::kernel)
		util_imported_link_libraries(algutil::Gecode algutil::Gecode::int)

		if(_GECODE_minimodel_LIBRARY)
			util_add_imported_library(algutil::Gecode::minimodel ${_GECODE_minimodel_LIBRARY} ${GECODE_INCLUDE_DIR})
			util_imported_link_libraries(algutil::Gecode::minimodel algutil::Gecode::int)
		endif()

		if(_GECODE_float_LIBRARY)
			util_add_imported_library(algutil::Gecode::float ${_GECODE_float_LIBRARY} ${GECODE_INCLUDE_DIR})
			util_imported_link_libraries(algutil::Gecode::float algutil::Gecode::int)
			util_imported_link_libraries(algutil::Gecode algutil::Gecode::float)

			if(_GECODE_minimodel_LIBRARY)
				util_imported_link_libraries(algutil::Gecode::minimodel algutil::Gecode::float)
			endif()
		endif()

		if(_GECODE_set_LIBRARY)
			util_add_imported_library(algutil::Gecode::set ${_GECODE_set_LIBRARY} ${GECODE_INCLUDE_DIR})
			util_imported_link_libraries(algutil::Gecode::set algutil::Gecode::int)
			util_imported_link_libraries(algutil::Gecode algutil::Gecode::set)

			if(_GECODE_minimodel_LIBRARY)
				util_imported_link_libraries(algutil::Gecode::minimodel algutil::Gecode::set)
			endif()
		endif()

		if(_GECODE_minimodel_LIBRARY)
			util_imported_link_libraries(algutil::Gecode algutil::Gecode::minimodel)
		endif()
	endif()

	if(_GECODE_search_LIBRARY)
		util_add_imported_library(algutil::Gecode::search ${_GECODE_search_LIBRARY} ${GECODE_INCLUDE_DIR})
		util_imported_link_libraries(algutil::Gecode::search algutil::Gecode::kernel)
		util_imported_link_libraries(algutil::Gecode algutil::Gecode::search)
	endif()

	if(_GECODE_search_LIBRARY AND _GECODE_int_LIBRARY AND _GECODE_gist_LIBRARY)
		util_add_imported_library(algutil::Gecode::gist ${_GECODE_gist_LIBRARY} ${GECODE_INCLUDE_DIR})
		util_imported_link_libraries(algutil::Gecode::gist algutil::Gecode::int)
		util_imported_link_libraries(algutil::Gecode::gist algutil::Gecode::search)
		util_imported_link_libraries(algutil::Gecode algutil::Gecode::gist)
	endif()

	if(_GECODE_driver_LIBRARY AND _GECODE_int_LIBRARY AND _GECODE_search_LIBRARY)
		util_add_imported_library(algutil::Gecode::driver ${_GECODE_driver_LIBRARY} ${GECODE_INCLUDE_DIR})
		util_imported_link_libraries(algutil::Gecode::driver algutil::Gecode::int)
		util_imported_link_libraries(algutil::Gecode::driver algutil::Gecode::search)
		if(_GECODE_minimodel_LIBRARY)
			util_imported_link_libraries(algutil::Gecode::driver algutil::Gecode::minimodel)
		endif()
		if(_GECODE_gist_LIBRARY)
			util_imported_link_libraries(algutil::Gecode::driver algutil::Gecode::gist)
		endif()
		util_imported_link_libraries(algutil::Gecode algutil::Gecode::driver)
	endif()

	if(_GECODE_flatzinc_LIBRARY)
		util_add_imported_library(algutil::Gecode::flatzinc ${_GECODE_flatzinc_LIBRARY} ${GECODE_INCLUDE_DIR})
		foreach(_LIB IN LISTS _GECODE_LIBS)
			if(NOT "${_LIB}" STREQUAL "" AND NOT "${_LIB}" STREQUAL "flatzinc")
				if(_GECODE_${_LIB}_LIBRARY)
					util_imported_link_libraries(algutil::Gecode::flatzinc "algutil::Gecode::${_LIB}")
				endif()
			endif()
		endforeach()
		util_imported_link_libraries(algutil::Gecode algutil::Gecode::flatzinc)
	endif()
endif(GECODE_FOUND)
endif(NOT TARGET algutil::Gecode)

