using SchrodingerSolver
using Base.Iterators
using GLMakie

savepath = "./scripts/stencil_plots/plots/"
basefilename = "stencil3d_"

for (i,j,k) in product(1:3,2:3,3:3)
    A = AssemblySymmetricOffset(AllZeroOffset, (1:i,1:j,1:k))
    f=Figure(size=(800,800))

    ax = Axis3(f[1,1])
    hidedecorations!(ax)
    hidespines!(ax)
    stencil!(ax,A)
    xlims!(ax,-2,2)
    ylims!(ax,-2,2)
    zlims!(ax,-2,2)
    save(savepath*basefilename*string(i)*"_"*string(j)*"_"*string(k)*".png",f)
end