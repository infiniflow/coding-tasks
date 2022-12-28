file(GLOB_RECURSE
        wavelet_src
        CONFIGURE_DEPENDS
        bit_array.cpp
        wavelet_tree.cpp
        wt_main.cpp
        )
add_executable(main ${wavelet_src})
target_link_libraries(main ${CMAKE_THREAD_LIBS_INIT})

