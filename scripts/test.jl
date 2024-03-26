using SchrodingerSolver

begin
    for storetype in [Int32, Int64]
        for floattype in [Float32, Float64]
            for ndims in 1:5
                begin
                    for h in range(floattype(0.2), floattype(1), 4)
                        ranges = ntuple(i -> range(floattype(3 + i), floattype(5 + i);
                                                   step=h),
                                        ndims)

                        grid = PeriodicGrid(storetype, floattype, floattype(0.01), ranges)

                        size(grid) == size.(ranges, 1)
                        length(grid) == prod(size.(ranges, 1))

                        begin
                            sz = size.(ranges, 1)
                            offset = CartesianIndex(size(grid))

                            println("Ranges: ", ranges)
                            println("Size: ", size(grid))

                            for idx in CartesianIndices(grid)
                                println(" ", getindex.(ranges, Tuple(idx)), " ", idx)

                                println(" ", getindex(grid, idx))
                            end

                            println("KS: ",
                                    [getindex.(ranges, Tuple(idx))
                                     for idx in CartesianIndices(grid)])

                            println("KV: ",
                                    [getindex(grid,
                                              idx + offset)
                                     for idx in CartesianIndices(grid)])

                            println("RS: ",
                                    [getindex(grid,
                                              idx + offset) ==
                                     getindex.(ranges, Tuple(idx))
                                     for idx in CartesianIndices(grid)])

                            all([getindex(grid, idx) ==
                                 getindex.(ranges, Tuple(idx))
                                 for idx in CartesianIndices(grid)])

                            all([getindex(grid,
                                          idx + offset) ==
                                 getindex.(ranges, Tuple(idx))
                                 for idx in CartesianIndices(grid)])

                            all([getindex(grid,
                                          idx - offset) ==
                                 getindex.(ranges, Tuple(idx))
                                 for idx in CartesianIndices(grid)])
                        end
                    end
                end
            end
        end
    end
end
