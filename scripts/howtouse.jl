using SchrodingerSolver

#Always define functions in a vectorial way
function f1(x)
    @views  2*x[:,1]-3*x[:,2]
end

function f2(x)
    @views 3*x[:,1]+4*x[:,2]
end

function psi1(x)
    @views 3*x[:,1]+5*x[:,2]
end

function psi2(x)
    @views 2*x[:,1]+4*x[:,2]
end

C1 = SchrodingerPDEComponent(0.5,f1,psi1)
C2 = SchrodingerPDEComponent(0.5,f2,psi1)
C3 = SchrodingerPDEComponent(0.5,f1,psi2)

function F(x)
    @views 2*x[:,1]+3*x[:,2]+4*x[:,3]
end

PDE = SchrodingerPDENonPolynomic(((1.0,2.0),(1.0,2.0)),(C1,C2,C3),F,1.0)

Mesh = PeriodicGrid(CPUBackendF64,PDE,1.0,(100,100))

Params = DefaultSolver(CPUBackendF64,2)

Method = PaulMethod1(CPUBackendF64,PDE,Mesh,Params)

#Problem = SchrodingerProblem(PDE,DefaultSolver(2),Mesh)

#step!(Problem)

