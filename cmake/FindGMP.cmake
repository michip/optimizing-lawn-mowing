include(FindPackageHandleStandardArgs)

if (NOT TARGET algutil::GMP AND NOT TARGET algutil::MPFR) # guard against multiple searches

    set(GMP_ROOT "" CACHE STRING "A path where the GNU multiprecision library (GMP) can be found; if unset, it is searched for in the system default directories.")
    set(MPFR_ROOT "" CACHE STRING "A path where the MPFR library can be found; if unset, it is searched for using GMP_ROOT and the system default directories.")
    mark_as_advanced(GMP_ROOT MPFR_ROOT)

    # unfortunately, there is no include directory for GMP except for the toplevel directory...
    find_path(
            GMP_INCLUDE_DIR
            NAMES
            gmp.h
            gmpxx.h
            HINTS
                /usr/local/Cellar/gmp/6.2.1/include
                ${GMP_ROOT}/include
                ${GMP_ROOT}
            PATH_SUFFIXES
            GMP
            gmp
            Gmp
    )
    mark_as_advanced(GMP_INCLUDE_DIR)

    find_library(
            GMP_LIBRARY
            NAMES
            "gmp"
            "gmp-9"
            "gmp-10"
            "gmp-11"
            "gmp-12"
            "libgmp.lib"
            "libgmp-9.lib"
            "libgmp-10.lib"
            "libgmp-11.lib"
            "libgmp-12.lib"
            HINTS
            /usr/local/Cellar/gmp/6.2.1/lib
            ${GMP_ROOT}/lib
            ${GMP_ROOT}
            PATH_SUFFIXES
            GMP
            gmp
            Gmp
    )
    mark_as_advanced(GMP_LIBRARY)

    find_library(
            MPFR_LIBRARY
            NAMES
            "mpfr"
            "mpfr-3"
            "mpfr-4"
            "mpfr-5"
            "mpfr-6"
            "libmpfr.lib"
            "libmpfr-3.lib"
            "libmpfr-4.lib"
            "libmpfr-5.lib"
            "libmpfr-6.lib"
            HINTS
            /usr/local/Cellar/mpfr/4.1.0/lib
            ${MPFR_ROOT}/lib
            ${MPFR_ROOT}
            ${GMP_ROOT}/lib
            ${GMP_ROOT}
            PATH_SUFFIXES
            MPFR
            mpfr
            Mpfr
    )
    mark_as_advanced(MPFR_LIBRARY)

    find_path(
            MPFR_INCLUDE_DIR
            NAMES
            mpfr.h
            HINTS
            /usr/local/Cellar/mpfr/4.1.0/include
            ${MPFR_ROOT}/include
            ${MPFR_ROOT}
            ${GMP_ROOT}/include
            ${GMP_ROOT}
            PATH_SUFFIXES
            MPFR
            mpfr
            Mpfr
    )
    mark_as_advanced(MPFR_INCLUDE_DIR)

    if (GMP_LIBRARY AND GMP_INCLUDE_DIR AND NOT TARGET algutil::GMP)
        util_add_imported_library(algutil::GMP ${GMP_LIBRARY} ${GMP_INCLUDE_DIR})

        if (MPFR_INCLUDE_DIR AND MPFR_LIBRARY AND NOT TARGET algutil::MPFR)
            util_add_imported_library(algutil::MPFR ${MPFR_LIBRARY} ${MPFR_INCLUDE_DIR})
            util_imported_link_libraries(algutil::MPFR algutil::GMP)
        endif ()
    endif ()

    find_package_handle_standard_args(GMP DEFAULT_MSG GMP_LIBRARY GMP_INCLUDE_DIR MPFR_LIBRARY MPFR_INCLUDE_DIR)

endif ()
