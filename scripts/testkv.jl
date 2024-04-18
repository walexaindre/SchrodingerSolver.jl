using SchrodingerSolver

mesh2d = PeriodicAbstractMesh(Int64, (10, 8))

A = get_A_format_IJV(Float64,mesh2d,(:ord2,:ord4))


A = sparse(A...)
drop(A,mesh2d,0.0001)

