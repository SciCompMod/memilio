if (TARGET epidemiology_io)
    add_executable(paper_202011 paper_202011.cpp)
    target_link_libraries(paper_202011 PRIVATE epidemiology epidemiology_io Boost::filesystem)
    target_compile_options(paper_202011 PRIVATE
        $<$<CXX_COMPILER_ID:GNU>: -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated-copy -Wno-unknown-warning-option -Wno-unknown-warning-option>
        $<$<CXX_COMPILER_ID:CLANG>: -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated -Wno-unknown-warning-option -Wno-unknown-warning-option>
        $<$<CXX_COMPILER_ID:MSVC>: /W4 /WX>)
endif (TARGET epidemiology_io)