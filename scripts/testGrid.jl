using SchrodingerSolver

import SchrodingerSolver: PeriodicGrid, get_D_format_IJV, PeriodicAbstractMesh, PeriodicGrid

pgrid = PeriodicGrid(Int, Float64, 3.5, 1.0:0.5:60.5,1:0.5:70.5,1:0.5:70.5)
I,J,V = get_D_format_IJV(Float64, pgrid,(:ord6,:ord8,:ord10)) 

