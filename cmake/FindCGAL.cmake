include(FindPackageHandleStandardArgs)

find_package(GMP)
find_package(Boost COMPONENTS system thread)

if(GMP_FOUND AND Boost_FOUND AND NOT TARGET algutil::CGAL)

    find_path(
            CGAL_INCLUDE_DIR
            NAMES
            CGAL/Arrangement_2.h
            HINTS
            /opt/tmp/include
            ${CGAL_ROOT}/include
            "C:/Program Files/CGAL/*/include"
            "C:/Program Files/CGAL/include"
            PATH_SUFFIXES
            CGAL
            Cgal
            cgal
    )

    function(util_cgal_construct_libstring VARNAME DEBUG RELEASE)
        if(DEBUG AND NOT RELEASE)
            set(${VARNAME} ${DEBUG} CACHE FILEPATH "Location of this CGAL library.")
        elseif(RELEASE AND NOT DEBUG)
            set(${VARNAME} ${RELEASE} CACHE FILEPATH "Location of this CGAL library.")
        elseif(RELEASE AND DEBUG)
            set(${VARNAME} "debug;${DEBUG};optimized;${RELEASE}" CACHE FILEPATH "Location of this CGAL library.")
        endif()
        if(NOT CGAL_FIND_QUIETLY)
            message(STATUS "${VARNAME}: ${${VARNAME}}")
        endif()
    endfunction(util_cgal_construct_libstring)

    if(MSVC AND CGAL_INCLUDE_DIR)
        set(_CGAL_INCLUDE_ROOT "${CGAL_INCLUDE_DIR}/..")

        file(GLOB _CGAL_WINDOWS_LIB_GLOB "${_CGAL_INCLUDE_ROOT}/lib/*.lib")
        foreach(F IN LISTS _CGAL_WINDOWS_LIB_GLOB)
            string(REGEX MATCH "/[^/]*-(vc[1-9][0-9][0-9])-[^/]*-([0-9\\.]+)\\.[Ll][Ii][Bb]$" _CGAL_WINDOWS_LIB_GLOB_RESULT ${F})
            if(_CGAL_WINDOWS_LIB_GLOB_RESULT)
                set(_CGAL_WINDOWS_LIB_TOOLSET ${CMAKE_MATCH_1})
                set(_CGAL_WINDOWS_LIB_VERSION ${CMAKE_MATCH_2})
                break()
            endif()
        endforeach()

        find_library(
                CGAL_LIBRARY_RELEASE
                NAMES
                "CGAL-${_CGAL_WINDOWS_LIB_TOOLSET}-mt-${_CGAL_WINDOWS_LIB_VERSION}.lib"
                HINTS
                "/opt/tmp/lib ${_CGAL_INCLUDE_ROOT}/lib"
        )

        find_library(
                CGAL_LIBRARY_DEBUG
                NAMES
                "CGAL-${_CGAL_WINDOWS_LIB_TOOLSET}-mt-gd-${_CGAL_WINDOWS_LIB_VERSION}.lib"
                HINTS
                /opt/tmp/lib
                ${_CGAL_INCLUDE_ROOT}/lib
        )

        mark_as_advanced(CGAL_LIBRARY_DEBUG CGAL_LIBRARY_RELEASE)
        util_cgal_construct_libstring(CGAL_LIBRARY "${CGAL_LIBRARY_DEBUG}" "${CGAL_LIBRARY_RELEASE}")

        find_library(
                CGAL_CORE_LIBRARY_RELEASE
                NAMES
                "CGAL_Core-${_CGAL_WINDOWS_LIB_TOOLSET}-mt-${_CGAL_WINDOWS_LIB_VERSION}.lib"
                HINTS
                /opt/tmp/lib
                ${_CGAL_INCLUDE_ROOT}/lib
        )

        find_library(
                CGAL_CORE_LIBRARY_DEBUG
                NAMES
                "CGAL_Core-${_CGAL_WINDOWS_LIB_TOOLSET}-mt-gd-${_CGAL_WINDOWS_LIB_VERSION}.lib"
                HINTS
                /opt/tmp/lib
                ${_CGAL_INCLUDE_ROOT}/lib
        )

        mark_as_advanced(CGAL_CORE_LIBRARY_DEBUG CGAL_CORE_LIBRARY_RELEASE)
        util_cgal_construct_libstring(CGAL_CORE_LIBRARY "${CGAL_CORE_LIBRARY_DEBUG}" "${CGAL_CORE_LIBRARY_RELEASE}")
    endif()

    if(NOT CGAL_LIBRARY OR NOT CGAL_CORE_LIBRARY)
        find_library(
                CGAL_LIBRARY
                NAMES
                CGAL
                HINTS
                /opt/tmp/lib
                ${_CGAL_INCLUDE_ROOT}/lib
                ${_CGAL_INCLUDE_ROOT}/usr/lib
                ${_CGAL_INCLUDE_ROOT}/usr/local/lib
                ${CGAL_ROOT}/lib
                ${CGAL_ROOT}/usr/lib
                ${CGAL_ROOT}/usr/local/lib
                $ENV{ProgramFiles}/CGAL/*/lib
                $ENV{SystemDrive}/CGAL/*/lib
                PATH_SUFFIXES
                CGAL
                cgal
                Cgal
        )

        find_library(
                CGAL_CORE_LIBRARY
                NAMES
                CGAL_Core
                HINTS
                /opt/tmp/lib
                ${_CGAL_INCLUDE_ROOT}/lib
                ${CGAL_ROOT}/lib
                $ENV{ProgramFiles}/CGAL/*/lib
                $ENV{SystemDrive}/CGAL/*/lib
                PATH_SUFFIXES
                CGAL
                cgal
                Cgal
        )
    endif()

    mark_as_advanced(CGAL_INCLUDE_DIR CGAL_LIBRARY CGAL_CORE_LIBRARY)

    if(CGAL_INCLUDE_DIR AND CGAL_LIBRARY AND CGAL_CORE_LIBRARY)
        if(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            set(CGAL_COMPILE_FLAGS "-frounding-math")
            if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
                set(CGAL_COMPILE_FLAGS "${CGAL_COMPILE_FLAGS};-Wno-ignored-optimization-argument")
            endif()
        endif()

        if(MSVC)
            set(CGAL_COMPILE_FLAGS "/fp:strict")
        endif()

        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
            set(CGAL_COMPILE_FLAGS "-fp-model;strict")
        endif()

        mark_as_advanced(CGAL_COMPILE_FLAGS)

        set(_CGAL_POSSIBLY_IMAGEIO "CGAL_INCLUDE_DIR")
        if(CGAL_FIND_COMPONENTS AND "ImageIO" IN_LIST CGAL_FIND_COMPONENTS)
            set(_CGAL_POSSIBLY_IMAGEIO "CGAL_IMAGEIO_LIBRARY")
        endif()

        if(MSVC)
            find_library(
                    CGAL_IMAGEIO_LIBRARY_RELEASE
                    NAMES
                    "CGAL_ImageIO-${_CGAL_WINDOWS_LIB_TOOLSET}-mt-${_CGAL_WINDOWS_LIB_VERSION}.lib"
                    HINTS
                    "/opt/tmp/lib ${_CGAL_INCLUDE_ROOT}/lib"
            )

            find_library(
                    CGAL_IMAGEIO_LIBRARY_DEBUG
                    NAMES
                    "CGAL_ImageIO-${_CGAL_WINDOWS_LIB_TOOLSET}-mt-gd-${_CGAL_WINDOWS_LIB_VERSION}.lib"
                    HINTS
                    "/opt/tmp/lib ${_CGAL_INCLUDE_ROOT}/lib"
            )

            mark_as_advanced(CGAL_IMAGEIO_LIBRARY_DEBUG CGAL_IMAGEIO_LIBRARY_RELEASE)
            util_cgal_construct_libstring(CGAL_IMAGEIO_LIBRARY "${CGAL_IMAGEIO_LIBRARY_DEBUG}" "${CGAL_IMAGEIO_LIBRARY_RELEASE}")
        endif()

        if(NOT CGAL_IMAGEIO_LIBRARY)
            find_library(
                    CGAL_IMAGEIO_LIBRARY
                    NAMES
                    CGAL_IMAGEIO
                    CGAL_ImageIO
                    CGAL_ImageIo
                    CGAL_imageio
                    HINTS
                    /opt/tmp/lib
                    ${_CGAL_INCLUDE_ROOT}/lib
                    ${CGAL_ROOT}/lib
                    $ENV{ProgramFiles}/CGAL/*/lib
                    $ENV{SystemDrive}/CGAL/*/lib
                    PATH_SUFFIXES
                    CGAL
                    cgal
                    Cgal
            )
        endif()
        mark_as_advanced(CGAL_IMAGEIO_LIBRARY)

        if(CGAL_IMAGEIO_LIBRARY)
            set(CGAL_IMAGEIO_FOUND TRUE)
        endif()

        set(_CGAL_POSSIBLY_QT "CGAL_INCLUDE_DIR")
        if(CGAL_FIND_COMPONENTS AND "qt" IN_LIST CGAL_FIND_COMPONENTS)
            set(_CGAL_POSSIBLY_QT "CGAL_QT_LIBRARY")
        endif()

        if(MSVC)
            find_library(
                    CGAL_QT_LIBRARY_RELEASE
                    NAMES
                    "CGAL_Qt5-${_CGAL_WINDOWS_LIB_TOOLSET}-mt-${_CGAL_WINDOWS_LIB_VERSION}.lib"
                    HINTS
                    "/opt/tmp/lib ${_CGAL_INCLUDE_ROOT}/lib"
            )

            find_library(
                    CGAL_QT_LIBRARY_DEBUG
                    NAMES
                    "CGAL_Qt5-${_CGAL_WINDOWS_LIB_TOOLSET}-mt-gd-${_CGAL_WINDOWS_LIB_VERSION}.lib"
                    HINTS
                    "/opt/tmp/lib ${_CGAL_INCLUDE_ROOT}/lib"
            )

            mark_as_advanced(CGAL_QT_LIBRARY_DEBUG CGAL_QT_LIBRARY_RELEASE)
            util_cgal_construct_libstring(CGAL_QT_LIBRARY "${CGAL_QT_LIBRARY_DEBUG}" "${CGAL_QT_LIBRARY_RELEASE}")
        endif()

        if(NOT CGAL_QT_LIBRARY)
            find_library(
                    CGAL_QT_LIBRARY
                    NAMES
                    CGAL_Qt
                    CGAL_Qt4
                    CGAL_Qt5
                    HINTS
                    /opt/tmp/lib
                    ${_CGAL_INCLUDE_ROOT}/lib
                    ${CGAL_ROOT}/lib
                    $ENV{ProgramFiles}/CGAL/*/lib
                    $ENV{SystemDrive}/CGAL/*/lib
                    PATH_SUFFIXES
                    CGAL
                    cgal
                    Cgal
            )
        endif()
        mark_as_advanced(CGAL_QT_LIBRARY)

        if(CGAL_QT_LIBRARY)
            set(CGAL_QT_FOUND TRUE)
        endif()
    endif()
endif()


find_package_handle_standard_args(CGAL DEFAULT_MSG CGAL_LIBRARY CGAL_CORE_LIBRARY CGAL_INCLUDE_DIR ${_CGAL_POSSIBLY_IMAGEIO} ${_CGAL_POSSIBLY_QT})

if(CGAL_FOUND)
    util_add_imported_library(algutil::CGAL "${CGAL_LIBRARY}" "${CGAL_INCLUDE_DIR}")

    set_target_properties(algutil::CGAL PROPERTIES INTERFACE_COMPILE_OPTIONS "${CGAL_COMPILE_FLAGS}")
    util_imported_link_libraries(algutil::CGAL "algutil::MPFR")

    util_add_imported_library(algutil::CGAL::Core "${CGAL_CORE_LIBRARY}" "")
    util_imported_link_libraries(algutil::CGAL::Core algutil::CGAL)

    if(CGAL_QT_FOUND)
        util_add_imported_library(algutil::CGAL::Qt "${CGAL_QT_LIBRARY}" "")
        util_imported_link_libraries(algutil::CGAL::Qt algutil::CGAL::Core)
    endif()

    if(CGAL_IMAGEIO_FOUND)
        util_add_imported_library(algutil::CGAL::ImageIO "${CGAL_IMAGEIO_LIBRARY}" "")
        util_imported_link_libraries(algutil::CGAL::ImageIO algutil::CGAL::Core)
    endif()
endif()
