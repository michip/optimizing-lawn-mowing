# find a google ortools installation
include(FindPackageHandleStandardArgs)

if(NOT TARGET algutil::ortools)
	find_path(GFLAGS_INCLUDE_DIR
		NAMES
			gflags/gflags.h
		HINTS
			${ORTOOLS_ROOT}/include
			/usr/local/include
	)

	find_path(GLOG_INCLUDE_DIR
		NAMES
			glog/logging.h
		HINTS
			${ORTOOLS_ROOT}/include
			/usr/local/include
	)

	find_path(ABSL_INCLUDE_DIR
		NAMES
			absl/base/config.h
		HINTS
			${ORTOOLS_ROOT}/include
			${ABSL_ROOT}/include
			/usr/local/include
	)

	find_path(COIN_INCLUDE_DIR
		NAMES
			coin/CoinUtility.hpp
		HINTS
			${ORTOOLS_ROOT}/include
			${COIN_ROOT}/include
	)

	find_path(PROTOBUF_INCLUDE_DIR
		NAMES
			google/protobuf/any.h
		HINTS
			${ORTOOLS_ROOT}/include
			${PROTOBUF_ROOT}/include
	)

	find_path(ORTOOLS_INCLUDE_DIR
		NAMES
			ortools/base/hash.h
		HINTS
			${ORTOOLS_ROOT}/include
	)

	set(ORTOOLS_LIBRARIES
		absl_container
		absl_synchronization
		absl_throw_delegate
		absl_time
		absl_strings
		absl_bad_any_cast
		absl_failure_signal_handler
		absl_hash
		absl_examine_stack
		absl_stacktrace
		absl_variant
		absl_optional
		absl_symbolize
		absl_int128
		absl_dynamic_annotations
		absl_base
		absl_bad_optional_access
		absl_stack_consumption
		absl_spinlock_wait
		protobuf
		glog
		gflags
		CbcSolver
		Cbc
		OsiCbc
		Cgl
		Clp
		ClpSolver
		OsiClp
		Osi
		CoinUtils
		ortools
	)

	set(_LIB_VARS "")
	set(_ORTOOLS_LIBS "")
	set(_ORTOOLS_LIBS_FOUND TRUE)
	foreach(_LIB IN LISTS ORTOOLS_LIBRARIES)
		if(NOT "${_LIB}" STREQUAL "")
			set(_LIB_VARS "${_LIB_VARS};_${_LIB}_LIBRARY")
			find_library(_${_LIB}_LIBRARY
				NAMES
					${_LIB}
				HINTS
					${ORTOOLS_ROOT}/lib
					${PROTOBUF_ROOT}/lib
					${ABSL_ROOT}/lib
					${COIN_ROOT}/lib
					/usr/local/lib
			)

			if(NOT _${_LIB}_LIBRARY)
				set(_ORTOOLS_LIBS_FOUND FALSE)
				message(WARNING "Could not find library ${_LIB}!")
			else()
				set(_ORTOOLS_LIBS ${_${_LIB}_LIBRARY})
			endif()
		endif()
	endforeach()

	if(NOT ${_ORTOOLS_LIBS_FOUND})
		set(_ORTOOLS_LIBS "")
	endif()

	find_package_handle_standard_args(ORTOOLS DEFAULT_MSG GFLAGS_INCLUDE_DIR GLOG_INCLUDE_DIR ABSL_INCLUDE_DIR COIN_INCLUDE_DIR PROTOBUF_INCLUDE_DIR ORTOOLS_INCLUDE_DIR _ORTOOLS_LIBS)
	if(ORTOOLS_FOUND)
		set(_ORTOOLS_IMPORTED_TARGETS "")

		foreach(_LIB IN LISTS ORTOOLS_LIBRARIES)
			if(NOT ${_LIB} STREQUAL "ortools" AND NOT ${_LIB} STREQUAL "")
				util_add_imported_library(algutil::ortools_deps::${_LIB} "${_${_LIB}_LIBRARY}" "")
				set(_ORTOOLS_IMPORTED_TARGETS "${_ORTOOLS_IMPORTED_TARGETS};algutil::ortools_deps::${_LIB}")
			endif()
		endforeach()
		util_add_imported_library(algutil::ortools "${_ortools_LIBRARY}" "${GFLAGS_INCLUDE_DIR};${GLOG_INCLUDE_DIR};${ABSL_INCLUDE_DIR};${COIN_INCLUDE_DIR};${PROTOBUF_INCLUDE_DIR};${ORTOOLS_INCLUDE_DIR}")
		util_imported_link_libraries(algutil::ortools "${_ORTOOLS_IMPORTED_TARGETS}")
	endif()
endif()

