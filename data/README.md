# MEmilio Data #

This directory contains data used for real-world epidemiological simulations. Additional datasets (e.g., daily reported COVID‑19 data for Germany) can be added using the MEmilio Epidata Package [MEmilio Python library](../pycode/memilio-epidata/README.rst).

## Directory Structure
The data is organized by region. Within a region folder (e.g., Germany), you will find the following subdirectories:

- **contacts**: Contact matrices, i.e., the number of daily contacts between persons of different age groups.
- **mobility**: Number of daily commuters and other travelers moving between regional entities.
- **pydata**: Data generated using the `Epidata` package from the MEmilio library.
