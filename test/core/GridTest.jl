@testset "GridTest.jl" begin
    for storetype in [Int32, Int64]
        for floattype in [Float32, Float64]
            for ndims in 1:5
                @testset "PeriodicGrid{$ndims} {$floattype} - {$storetype}" begin
                    for h in range(floattype(0.2), floattype(1), 4)
                        ranges = ntuple(i -> range(floattype(3 + i), floattype(5 + i);
                                                   step=h),
                                        ndims)

                        grid = PeriodicGrid(storetype, floattype, floattype(0.01), ranges)

                        @test size(grid) == size.(ranges, 1)
                        @test length(grid) == prod(size.(ranges, 1))

                        @testset "Test getindex in Periodic boundaries" begin
                            sz = size.(ranges, 1)
                            offset = CartesianIndex(size(grid))
                            grid_indexes = CartesianIndices(grid)

                            @test all([getindex(grid, idx) ==
                                       getindex.(ranges, Tuple(idx))
                                       for idx in grid_indexes])

                            @test all([getindex(grid,
                                                idx + offset) ==
                                       getindex.(ranges, Tuple(idx))
                                       for idx in grid_indexes])

                            @test all([getindex(grid,
                                                idx - offset) ==
                                       getindex.(ranges, Tuple(idx))
                                       for idx in grid_indexes])
                        end
                    end
                end
            end
        end
    end
end
