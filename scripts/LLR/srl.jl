using GLMakie
using GeometryBasics: LineFace, GLIndex, Tesselation

fm(x) = (x[1]^2 + x[2] .^ 2 + 25 * (sin(x[1]) .^ 2 + sin(x[2]) .^ 2));
bounds = repeat([0 - 10.0, 10.0], 1, 2)
dens = 200
dX = (bounds[2, 1] - bounds[1, 1]) / dens
dY = (bounds[2, 2] - bounds[1, 2]) / dens
AX = bounds[1, 1]:dX:bounds[2, 1]
AY = bounds[1, 2]:dY:bounds[2, 2]
Z = fm.(tuple.(AX, AY'))

points = vec(Point3f.(AX, AY', Z))
# Connect the vetices with faces, as one would use for a 2D Rectangle
# grid with M,N grid points
faces = decompose(LineFace{GLIndex}, Tesselation(Rect2(0, 0, 1, 1), size(Z)))
grid_points = connect(points, faces)
colors = reinterpret(Float32, connect(Point{1, Float32}.(vec(Z)), faces))

linesegments(grid_points; color=collect(colors), colormap=:viridis, axis=(; type=Axis3), transparency=true, linewidth=0.5)