using SchrodingerSolver
using Base.Iterators
using GLMakie

savepath = "./scripts/stencil_plots/plots/"
basefilename = "stencil2d_"

for (i,j) in product(1:3,2:3)
    A = AssemblySymmetricOffset(AllZeroOffset, (1:i,1:j))
    f=Figure(size=(800,800))

    ax = Axis(f[1,1])
    hidedecorations!(ax)
    hidespines!(ax)
    stencil!(ax,A)
    xlims!(ax,-4,4)
    ylims!(ax,-4,4)
    save(savepath*basefilename*string(i)*"_"*string(j)*".png",f)
end