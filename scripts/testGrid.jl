using SchrodingerSolver

import SchrodingerSolver: PeriodicGrid, get_Δ_format_IJV, PeriodicAbstractMesh, PeriodicGrid

pgrid = PeriodicGrid(Int, Float64, 3.5, 1.0:0.5:3.0, 1.0:0.5:3.0)
get_Δ_format_IJV(Float64, pgrid,(:ord2,:ord4)) 

