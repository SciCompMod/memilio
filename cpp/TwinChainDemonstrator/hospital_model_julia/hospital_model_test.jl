using StatsBase

Base.@ccallable function printint(n::Cint)::Cvoid
    println("Your integer is $n")
    return nothing
end

Base.@ccallable function mean_of_two(x::Cint, y::Cint)::Cfloat
    return mean([x,y])
end