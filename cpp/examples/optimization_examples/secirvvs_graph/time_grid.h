#pragma once

#include <vector>
#include <cstddef>

template <typename FP>
std::vector<FP> make_time_grid(FP t0, FP tmax, size_t num_intervals)
{
    // ------------------------------------------------------------------- //
    // Create a uniform distribution of points in the specified timeframe. //
    // ------------------------------------------------------------------ //
    std::vector<FP> grid(num_intervals + 1);
    FP grid_spacing = (tmax - t0) / num_intervals;
    for (size_t i = 0; i <= num_intervals; i++) {
        grid[i] = t0 + i * grid_spacing;
    }
    return grid;
}