using SchrodingerSolver
using BenchmarkTools
using GLMakie

SUITE = BenchmarkGroup()
SUITE["rand"] = @benchmarkable rand(10)

# Write your benchmarks here.
tune!(SUITE)
run(SUITE)