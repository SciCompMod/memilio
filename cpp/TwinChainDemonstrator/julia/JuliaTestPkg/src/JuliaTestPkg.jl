module JuliaTestPkg
using StatsBase

export printint
export mean_of_two

# Julia mirror of the C struct in JuliaTestPkg.h:
# struct MyStruct { int a; double b; };
struct MyStruct
    a::Cint
    b::Cdouble
end

Base.@ccallable function printint(n::Cint)::Cvoid
    println("Your integer is $n")
    flush(stdout)
    return nothing
end

Base.@ccallable function mean_of_two(x::Cint, y::Cint)::Cfloat
    return mean([x,y])
end

Base.@ccallable function struct_test(s::MyStruct)::Cint
    println("My integer is $(s.a)")
    flush(stdout)
    return s.a
end

end # module JuliaTestPkg
