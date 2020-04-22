# This script installs package dependencies.

using Pkg

dependencies = [
    "Distributions",
    "Dates",
    "JLD",
    "HDF5",
    "DataFrames",
    "CSV",
    "Combinatorics",
    "Printf",
    "Distributed",
    "Random",
    "DelimitedFiles",
    "HypothesisTests",
]

Pkg.add(dependencies)