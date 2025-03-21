#=
Script to creat quick and dirty plots of a CSV file using julia's Makie visualization library. 
The input csv file is expected to have a column named `Time`. 
The standard julia environement is expected to have 
    - CSV
    - DataFrames
    - ColorSchemes
    - CairoMakie
installed. 
The plot is saved as a png file using the same filename as the input CSV.
=#


using CSV, DataFrames
using ColorSchemes
using CairoMakie

filename =ARGS[1]
data = CSV.read(filename, DataFrame)

fig = Figure()
ax = Axis(fig[1,1], title = filename, xlabel = "Time [d]")
for name in names(data)
    if name != "Time"
        lines!(data.Time, data[:, name], label = name, colormap=ColorSchemes.batlowKS)
    end
end
Legend(fig[1,2], ax)
save("$(filename[1:end-3])png", fig)
