using CSV, DataFrames, Plots

output_directory = "/home/nina/PycharmProjects/memilio/analysis/systems_files"
df = CSV.read(joinpath(output_directory, "concentration_curves.csv"), DataFrame, delim=",")


location = "MW043"
df_sub = df[(df.location.==location),:]
t = df_sub.time 
c = df_sub[!,"COV19"] 
p1 = plot(t/60, c, xlabel="t (h)", label=missing, ylabel="RNA concentration", thickness_scaling = 1.8)
savefig(p1, joinpath("analysis/plots/", "concentration_$location.png"))