# --------------------------------------------------------------------------- #
# PBeFS.cmake                                                                 #
#                                                                             #
# Fetches PBeFS (the library of parallel best-first searches whose            #
# branch-and-bound may explore the tree) and ADS (on which its parallel       #
# searches run) into the build tree, as bxs_fetch_pbefs(), after which        #
# pbefs_SOURCE_DIR and ads_SOURCE_DIR are the directories whose include/ is   #
# to add to the include paths. Only their headers are used, so the sources    #
# are fetched, not built: PBeFS at the commit BXS_PBEFS_COMMIT, ADS at the    #
# commit that this PBeFS pins. FETCHCONTENT_SOURCE_DIR_PBEFS and              #
# FETCHCONTENT_SOURCE_DIR_ADS give a source tree, or an installation, in      #
# place of the fetch, as a release tarball does [see the umbrella]:           #
#     bxs_fetch_pbefs()                                                       #
#     target_include_directories(<target> PRIVATE                             #
#                                ${pbefs_SOURCE_DIR}/include                  #
#                                ${ads_SOURCE_DIR}/include)                   #
# --------------------------------------------------------------------------- #

include(FetchContent)

set(BXS_PBEFS_COMMIT 3fc57057e54ca77edfb50d7d2dbe1a876f5ce12c
    CACHE STRING "The commit of PBeFS to download")

macro(bxs_fetch_pbefs)
    FetchContent_Declare(
            PBeFS
            GIT_REPOSITORY https://github.com/dmeoli/PBeFS.git
            GIT_TAG ${BXS_PBEFS_COMMIT}
            # the sources are fetched, not built
            SOURCE_SUBDIR none)
    FetchContent_MakeAvailable(PBeFS)

    # the commit of ADS is the one that the CMakeLists.txt of PBeFS pins,
    # which an installation of PBeFS does not have: then ADS has to be given
    set(_bxs_ads_commit "")
    if (EXISTS "${pbefs_SOURCE_DIR}/CMakeLists.txt")
        file(STRINGS "${pbefs_SOURCE_DIR}/CMakeLists.txt" _bxs_ads_pin
             REGEX "^[a-z_]+\\(PBEFS_ADS_COMMIT +[0-9a-f]+$")
        string(REGEX MATCH "[0-9a-f]+$" _bxs_ads_commit "${_bxs_ads_pin}")
    endif ()
    if (NOT _bxs_ads_commit AND NOT FETCHCONTENT_SOURCE_DIR_ADS)
        message(FATAL_ERROR "PBeFS in ${pbefs_SOURCE_DIR} pins no commit of "
                "ADS: give ADS by FETCHCONTENT_SOURCE_DIR_ADS")
    endif ()
    FetchContent_Declare(
            ADS
            GIT_REPOSITORY https://github.com/dmeoli/ADS.git
            GIT_TAG ${_bxs_ads_commit}
            SOURCE_SUBDIR none)
    FetchContent_MakeAvailable(ADS)
    unset(_bxs_ads_pin)
    unset(_bxs_ads_commit)
endmacro()
