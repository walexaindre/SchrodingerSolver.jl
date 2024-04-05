using SchrodingerSolver

import SchrodingerSolver: PeriodicGrid, get_Δ_format_IJV, PeriodicAbstractMesh, PeriodicGrid

pgrid = PeriodicGrid(Int, Float64, 3.5, 1.0:0.5:6.5)
I,J,V = get_Δ_format_IJV(Float64, pgrid,(:ord6,)) 

