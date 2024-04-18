using BenchmarkTools

struct Offset{M}
    offsets::NTuple{M, Int}
end

Base.length(o::Offset{M}) where{M}  = M

struct SymmetricOffsets{N}
    offsets::NTuple{N,Offset}
end

Base.length(o::SymmetricOffsets) = sum(length,o.offsets)

function apply_offsets(val::NTuple{N,Int},offsets::SymmetricOffsets{N}) where N

    res = fill(val,length(offsets))
    
    idx = 1
    dim = 1
    for offset in offsets.offsets
        for oidx in offset.offsets
            tmp = val[dim] + oidx
            res[idx] = Base.setindex(res[idx],tmp,dim)
            idx += 1
        end
        dim += 1
    end
    res
end


struct SymmetricOffsets2{M1,M2,M3}
    offsets::Tuple{Offset{M1},Offset{M2},Offset{M3}}
end


Base.length(o::SymmetricOffsets2{M1,M2,M3}) where {M1,M2,M3} = M1+M2+M3

function apply_offsets(val::NTuple{N,Int},offsets::SymmetricOffsets2{M1,M2,M3}) where {N,M1,M2,M3}

    res = fill(val,length(offsets))
    
    idx = 1
    for dim in 1:3
        for oidx in offsets.offsets[dim].offsets
            tmp = val[dim] + oidx
            res[idx] = Base.setindex(res[idx],tmp,dim)
            idx += 1
        end
    end

    res
end

struct SymmetricOffsets3{M1,M2,M3}
    o1::Offset{M1}
    o2::Offset{M2}
    o3::Offset{M3}
end

Base.length(o::SymmetricOffsets3{M1,M2,M3}) where {M1,M2,M3} = M1+M2+M3


function apply_offsets(val::NTuple{N,Int},offsets::SymmetricOffsets3{M1,M2,M3}) where {N,M1,M2,M3}

    res = fill(val,length(offsets))
    
    idx = 1

    for oidx in offsets.o1.offsets
        tmp = val[1] + oidx
        res[idx] = Base.setindex(res[idx],tmp,1)
        idx += 1
    end

    for oidx in offsets.o2.offsets
        tmp = val[2] + oidx
        res[idx] = Base.setindex(res[idx],tmp,2)
        idx += 1
    end

    for oidx in offsets.o3.offsets
        tmp = val[3] + oidx
        res[idx] = Base.setindex(res[idx],tmp,3)
        idx += 1
    end

    res
end




offset1 = Offset((1,2,3))
offset2 = Offset((4,5,6,7))
offset3 = Offset((8,9,10,11,12))

symmetric_offsets = SymmetricOffsets((offset1,offset2,offset3))
symmetric_offsets2 = SymmetricOffsets2((offset1,offset2,offset3))
symmetric_offsets3 = SymmetricOffsets3(offset1,offset2,offset3)


@btime apply_offsets($(1,2,3),$symmetric_offsets)
@btime apply_offsets($(1,2,3),$symmetric_offsets2)
@btime apply_offsets($(1,2,3),$symmetric_offsets3)

@code_warntype apply_offsets((1,2,3),symmetric_offsets)
@code_warntype apply_offsets((1,2,3),symmetric_offsets2)
@code_warntype apply_offsets((1,2,3),symmetric_offsets3)