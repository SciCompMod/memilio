module JuliaTestPkg
using StatsBase
using EnumX

export printint
export mean_of_two

# Julia mirror of the C struct in JuliaTestPkg.h:
# struct MyStruct { int a; double b; };
struct MyStruct
    a::Cint
    b::Cdouble
end

@enumx InfectionState::UInt32 begin
    Susceptible = 0
    Exposed
    InfectedNoSymptoms
    InfectedSymptoms
    InfectedSevere
    InfectedCritical
    Recovered
    Dead
    Count
end

@enumx ActivityType::UInt32 begin
    Home             = 0
    School              # 1
    Work                # 2
    Recreation          # 3
    BasicsShop          # 4
    Hospital            # 5
    ICU                 # 6
    Cemetery            # 7
    PublicTransport     # 8
    Count               # 9
    Invalid             # 10
end

struct InfectionTuple
    time_point::Cint                # time point of aquiring the infection state in seconds since the start of the simulation
    infection_state::InfectionState.T
end

struct PersonData
    id::Cint
    age_group::Cint
    activity_type::ActivityType.T
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

Base.@ccallable function hospital_model(p::PersonData,
    history::Ptr{InfectionTuple},
    history_len::Csize_t,
    time_of_day::Cint,
)::Bool
    println("Person id: $(p.id)")
    flush(stdout)
    if time_of_day >= 18
        return true
    else
        return false
    end
end

end # module JuliaTestPkg
