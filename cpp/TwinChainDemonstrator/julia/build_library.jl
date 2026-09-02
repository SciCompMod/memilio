using PackageCompiler

create_library(
    joinpath(@__DIR__, "JuliaTestPkg"),
    joinpath(@__DIR__, "..", "julia_lib");
    lib_name = "juliatestpkg",
    precompile_execution_file = joinpath(@__DIR__, "precompile.jl"),
    header_files = [joinpath(@__DIR__, "JuliaTestPkg.h")],
    incremental = false,
    filter_stdlibs = false,
    force = true,
)
