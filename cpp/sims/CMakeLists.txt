if (TARGET epidemiology_io)
    add_executable(paper_202011 paper_202011.cpp)
    target_link_libraries(paper_202011 PRIVATE epidemiology epidemiology_io Boost::filesystem)
    target_compile_options(paper_202011 PRIVATE
        $<$<CXX_COMPILER_ID:GNU>: -Wno-unknown-warning -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated-copy>
        $<$<CXX_COMPILER_ID:CLANG>: -Wno-unknown-warning-option -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated>
        $<$<CXX_COMPILER_ID:MSVC>: /W4 /WX>)
endif (TARGET epidemiology_io)