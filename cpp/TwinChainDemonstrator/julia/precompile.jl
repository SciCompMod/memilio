using JuliaTestPkg
JuliaTestPkg.printint(Int32(42))
JuliaTestPkg.mean_of_two(Int32(1), Int32(2))
JuliaTestPkg.struct_test(JuliaTestPkg.MyStruct(10, 3.14))
JuliaTestPkg.hospital_model(JuliaTestPkg.PersonData(0, 0, JuliaTestPkg.ActivityType.Home), Ptr{JuliaTestPkg.InfectionTuple}(C_NULL), Csize_t(0), Cint(20))
