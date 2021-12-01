include(FindPackageHandleStandardArgs)

if(NOT TARGET algutil::ODBC) # guard against multiple searches
	if(NOT WIN32)
		set(ODBC_INCLUDE_DIR "" CACHE STRING "The directory where the headers (sql.h etc.) for the ODBC library can be found.")
		set(ODBC_LIBRARY_DIR "" CACHE STRING "The directory where the ODBC library can be found.")
		mark_as_advanced(ODBC_INCLUDE_DIR ODBC_LIBRARY_DIR)

		find_path(
			ODBC_INCLUDE_DIR_
			NAMES
				sql.h
				sqlext.h
			HINTS
				${ODBC_INCLUDE_DIR}
				${ODBC_INCLUDE_DIR}/include
				/opt/local/include
				/usr/local/include
			PATH_SUFFIXES
				odbc
				ODBC
				Odbc
		)

		find_library(
			ODBC_LIBRARY_
			NAMES
				odbc
				libodbc.a
				odbc32
				libodbc32.lib
				odbc64
				libodbc64.lib
			HINTS
				${ODBC_LIBRARY_DIR}
				${ODBC_LIBRARY_DIR}/lib
				${ODBC_LIBRARY_DIR}/bin
				/opt/local/lib
				/opt/local/bin
				/usr/local/lib
				/usr/local/bin
			PATH_SUFFIXES
				odbc
				ODBC
				Odbc
		)

		find_package_handle_standard_args(ODBC DEFAULT_MSG ODBC_LIBRARY_ ODBC_INCLUDE_DIR_)
		set(ODBC_INCLUDE_DIRS ${ODBC_INCLUDE_DIR_})
		set(ODBC_LIBRARIES ${ODBC_LIBRARY_})
		util_add_imported_library(algutil::ODBC ${ODBC_LIBRARY_} ${ODBC_INCLUDE_DIR_})
	else()
		util_add_imported_library(algutil::ODBC "odbc32.lib" "")
		set(ODBC_INCLUDE_DIRS "")
		set(ODBC_LIBRARIES "odbc32.lib")
		find_package_handle_standard_args(ODBC DEFAULT_MSG ODBC_LIBRARIES)
	endif()
endif()

