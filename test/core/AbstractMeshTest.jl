@testset "AbstractMesh.jl" begin
    for ndims in 1:5
        for type in [Int32, Int64]
            @testset "PeriodicAbstractMesh{$ndims} {$type}" begin
                for dims in product(ntuple(i -> type(3):type(5), ndims)...)
                    mesh = PeriodicAbstractMesh(type, dims)

                    @test size(mesh) == dims
                    @test length(mesh) == prod(dims)

                    @testset "Test getindex in Periodic boundaries" begin
                        @test all(getindex(mesh, idx) == LinearIndices(dims)[idx]
                                  for idx in CartesianIndices(mesh))
                        @test all(getindex(mesh, idx + CartesianIndex(dims)) ==
                                  LinearIndices(dims)[idx]
                                  for idx in CartesianIndices(mesh))
                        @test all(getindex(mesh, idx - CartesianIndex(dims)) ==
                                  LinearIndices(dims)[idx]
                                  for idx in CartesianIndices(mesh))
                    end
                end
            end
        end
    end
end