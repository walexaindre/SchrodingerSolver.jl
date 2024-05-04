#Compares: 29 Depth 8 
function swapsort(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10;lt=isless,by=identity)
    i1,i8 = min_max(i1,i8,lt,by)
    i2,i7 = min_max(i2,i7,lt,by)
    i3,i10 = min_max(i3,i10,lt,by)
    i4,i9 = min_max(i4,i9,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)

    i1,i4 = min_max(i1,i4,lt,by)
    i2,i5 = min_max(i2,i5,lt,by)
    i6,i9 = min_max(i6,i9,lt,by)
    i7,i10 = min_max(i7,i10,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i4,i7 = min_max(i4,i7,lt,by)
    i8,i10 = min_max(i8,i10,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)
    i6,i8 = min_max(i6,i8,lt,by)
    i9,i10 = min_max(i9,i10,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i3,i6 = min_max(i3,i6,lt,by)
    i5,i8 = min_max(i5,i8,lt,by)
    i7,i9 = min_max(i7,i9,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)

    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    return i1, i2, i3, i4, i5, i6, i7, i8, i9, i10
end

swapsort(x::NTuple{10,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10],lt=lt,by=by)


#Compares: 35 Depth 8 
function swapsort(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11;lt=isless,by=identity)
    i1,i11 = min_max(i1,i11,lt,by)
    i2,i9 = min_max(i2,i9,lt,by)
    i3,i10 = min_max(i3,i10,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)

    i1,i7 = min_max(i1,i7,lt,by)
    i2,i4 = min_max(i2,i4,lt,by)
    i5,i10 = min_max(i5,i10,lt,by)
    i6,i8 = min_max(i6,i8,lt,by)
    i9,i11 = min_max(i9,i11,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)
    i4,i7 = min_max(i4,i7,lt,by)
    i6,i9 = min_max(i6,i9,lt,by)
    i8,i11 = min_max(i8,i11,lt,by)

    i3,i6 = min_max(i3,i6,lt,by)
    i4,i10 = min_max(i4,i10,lt,by)
    i5,i9 = min_max(i5,i9,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i2,i5 = min_max(i2,i5,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i7,i9 = min_max(i7,i9,lt,by)
    i8,i10 = min_max(i8,i10,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)
    i6,i8 = min_max(i6,i8,lt,by)
    i10,i11 = min_max(i10,i11,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)
    i9,i10 = min_max(i9,i10,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)
    return i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11
end

swapsort(x::NTuple{11,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11],lt=lt,by=by)


#Compares: 39 Depth 9 
function swapsort(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12;lt=isless,by=identity)
    i1,i12 = min_max(i1,i12,lt,by)
    i2,i11 = min_max(i2,i11,lt,by)
    i3,i10 = min_max(i3,i10,lt,by)
    i4,i9 = min_max(i4,i9,lt,by)
    i5,i8 = min_max(i5,i8,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)

    i1,i6 = min_max(i1,i6,lt,by)
    i2,i4 = min_max(i2,i4,lt,by)
    i7,i12 = min_max(i7,i12,lt,by)
    i9,i11 = min_max(i9,i11,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i4,i8 = min_max(i4,i8,lt,by)
    i5,i9 = min_max(i5,i9,lt,by)
    i10,i12 = min_max(i10,i12,lt,by)

    i2,i5 = min_max(i2,i5,lt,by)
    i3,i6 = min_max(i3,i6,lt,by)
    i7,i10 = min_max(i7,i10,lt,by)
    i8,i11 = min_max(i8,i11,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)
    i4,i7 = min_max(i4,i7,lt,by)
    i6,i9 = min_max(i6,i9,lt,by)
    i8,i10 = min_max(i8,i10,lt,by)
    i11,i12 = min_max(i11,i12,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i5,i8 = min_max(i5,i8,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    i9,i11 = min_max(i9,i11,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i7,i9 = min_max(i7,i9,lt,by)
    i10,i11 = min_max(i10,i11,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)
    i9,i10 = min_max(i9,i10,lt,by)

    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)
    return i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12
end

swapsort(x::NTuple{12,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12],lt=lt,by=by)


#Compares: 45 Depth 10 
function swapsort(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13;lt=isless,by=identity)
    i1,i12 = min_max(i1,i12,lt,by)
    i2,i13 = min_max(i2,i13,lt,by)
    i3,i11 = min_max(i3,i11,lt,by)
    i4,i10 = min_max(i4,i10,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)
    i6,i9 = min_max(i6,i9,lt,by)

    i1,i4 = min_max(i1,i4,lt,by)
    i2,i6 = min_max(i2,i6,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)
    i7,i11 = min_max(i7,i11,lt,by)
    i9,i13 = min_max(i9,i13,lt,by)
    i10,i12 = min_max(i10,i12,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i4,i9 = min_max(i4,i9,lt,by)
    i6,i10 = min_max(i6,i10,lt,by)
    i8,i12 = min_max(i8,i12,lt,by)
    i11,i13 = min_max(i11,i13,lt,by)

    i5,i10 = min_max(i5,i10,lt,by)
    i7,i9 = min_max(i7,i9,lt,by)
    i8,i11 = min_max(i8,i11,lt,by)
    i12,i13 = min_max(i12,i13,lt,by)

    i2,i8 = min_max(i2,i8,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    i9,i11 = min_max(i9,i11,lt,by)
    i10,i12 = min_max(i10,i12,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i8 = min_max(i3,i8,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i9,i10 = min_max(i9,i10,lt,by)
    i11,i12 = min_max(i11,i12,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i5,i8 = min_max(i5,i8,lt,by)
    i7,i9 = min_max(i7,i9,lt,by)
    i10,i11 = min_max(i10,i11,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i3,i6 = min_max(i3,i6,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)
    i9,i10 = min_max(i9,i10,lt,by)

    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    return i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13
end

swapsort(x::NTuple{13,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13],lt=lt,by=by)


#Compares: 51 Depth 10 
function swapsort(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14;lt=isless,by=identity)
    i1,i14 = min_max(i1,i14,lt,by)
    i2,i13 = min_max(i2,i13,lt,by)
    i3,i12 = min_max(i3,i12,lt,by)
    i4,i11 = min_max(i4,i11,lt,by)
    i5,i10 = min_max(i5,i10,lt,by)
    i6,i9 = min_max(i6,i9,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)

    i1,i5 = min_max(i1,i5,lt,by)
    i2,i6 = min_max(i2,i6,lt,by)
    i4,i7 = min_max(i4,i7,lt,by)
    i8,i11 = min_max(i8,i11,lt,by)
    i9,i13 = min_max(i9,i13,lt,by)
    i10,i14 = min_max(i10,i14,lt,by)

    i1,i4 = min_max(i1,i4,lt,by)
    i3,i10 = min_max(i3,i10,lt,by)
    i5,i12 = min_max(i5,i12,lt,by)
    i11,i14 = min_max(i11,i14,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i5,i8 = min_max(i5,i8,lt,by)
    i7,i10 = min_max(i7,i10,lt,by)
    i9,i11 = min_max(i9,i11,lt,by)
    i12,i13 = min_max(i12,i13,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i7 = min_max(i3,i7,lt,by)
    i5,i9 = min_max(i5,i9,lt,by)
    i6,i10 = min_max(i6,i10,lt,by)
    i8,i12 = min_max(i8,i12,lt,by)
    i13,i14 = min_max(i13,i14,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)
    i6,i9 = min_max(i6,i9,lt,by)
    i10,i12 = min_max(i10,i12,lt,by)
    i11,i13 = min_max(i11,i13,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    i6,i8 = min_max(i6,i8,lt,by)
    i7,i9 = min_max(i7,i9,lt,by)
    i10,i11 = min_max(i10,i11,lt,by)
    i12,i13 = min_max(i12,i13,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)
    i8,i10 = min_max(i8,i10,lt,by)
    i11,i12 = min_max(i11,i12,lt,by)

    i5,i6 = min_max(i5,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)
    i9,i10 = min_max(i9,i10,lt,by)

    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)
    i10,i11 = min_max(i10,i11,lt,by)
    return i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14
end

swapsort(x::NTuple{14,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14],lt=lt,by=by)


#Compares: 56 Depth 10 
function swapsort(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15;lt=isless,by=identity)
    i1,i13 = min_max(i1,i13,lt,by)
    i2,i10 = min_max(i2,i10,lt,by)
    i3,i11 = min_max(i3,i11,lt,by)
    i5,i12 = min_max(i5,i12,lt,by)
    i6,i15 = min_max(i6,i15,lt,by)
    i7,i14 = min_max(i7,i14,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)

    i1,i5 = min_max(i1,i5,lt,by)
    i3,i7 = min_max(i3,i7,lt,by)
    i4,i10 = min_max(i4,i10,lt,by)
    i6,i8 = min_max(i6,i8,lt,by)
    i9,i15 = min_max(i9,i15,lt,by)
    i11,i14 = min_max(i11,i14,lt,by)
    i12,i13 = min_max(i12,i13,lt,by)

    i1,i6 = min_max(i1,i6,lt,by)
    i2,i9 = min_max(i2,i9,lt,by)
    i4,i8 = min_max(i4,i8,lt,by)
    i5,i11 = min_max(i5,i11,lt,by)
    i7,i12 = min_max(i7,i12,lt,by)
    i10,i14 = min_max(i10,i14,lt,by)
    i13,i15 = min_max(i13,i15,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i2,i5 = min_max(i2,i5,lt,by)
    i4,i7 = min_max(i4,i7,lt,by)
    i8,i12 = min_max(i8,i12,lt,by)
    i9,i11 = min_max(i9,i11,lt,by)
    i10,i13 = min_max(i10,i13,lt,by)
    i14,i15 = min_max(i14,i15,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i3,i6 = min_max(i3,i6,lt,by)
    i5,i8 = min_max(i5,i8,lt,by)
    i7,i9 = min_max(i7,i9,lt,by)
    i11,i13 = min_max(i11,i13,lt,by)
    i12,i14 = min_max(i12,i14,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i4 = min_max(i3,i4,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)
    i6,i10 = min_max(i6,i10,lt,by)
    i11,i12 = min_max(i11,i12,lt,by)
    i13,i14 = min_max(i13,i14,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i8,i10 = min_max(i8,i10,lt,by)
    i9,i11 = min_max(i9,i11,lt,by)
    i12,i13 = min_max(i12,i13,lt,by)

    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)
    i10,i11 = min_max(i10,i11,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)
    i9,i10 = min_max(i9,i10,lt,by)
    i11,i12 = min_max(i11,i12,lt,by)

    i6,i7 = min_max(i6,i7,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)
    return i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15
end

swapsort(x::NTuple{15,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15],lt=lt,by=by)


#Compares: 60 Depth 10 
function swapsort(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16;lt=isless,by=identity)
    i1,i16 = min_max(i1,i16,lt,by)
    i2,i15 = min_max(i2,i15,lt,by)
    i3,i14 = min_max(i3,i14,lt,by)
    i4,i13 = min_max(i4,i13,lt,by)
    i5,i12 = min_max(i5,i12,lt,by)
    i6,i11 = min_max(i6,i11,lt,by)
    i7,i10 = min_max(i7,i10,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)

    i1,i8 = min_max(i1,i8,lt,by)
    i2,i7 = min_max(i2,i7,lt,by)
    i3,i6 = min_max(i3,i6,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    i9,i16 = min_max(i9,i16,lt,by)
    i10,i15 = min_max(i10,i15,lt,by)
    i11,i14 = min_max(i11,i14,lt,by)
    i12,i13 = min_max(i12,i13,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i2,i4 = min_max(i2,i4,lt,by)
    i5,i10 = min_max(i5,i10,lt,by)
    i6,i8 = min_max(i6,i8,lt,by)
    i7,i12 = min_max(i7,i12,lt,by)
    i9,i11 = min_max(i9,i11,lt,by)
    i13,i15 = min_max(i13,i15,lt,by)
    i14,i16 = min_max(i14,i16,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i7 = min_max(i3,i7,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i5,i9 = min_max(i5,i9,lt,by)
    i8,i12 = min_max(i8,i12,lt,by)
    i10,i14 = min_max(i10,i14,lt,by)
    i11,i13 = min_max(i11,i13,lt,by)
    i15,i16 = min_max(i15,i16,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)
    i6,i11 = min_max(i6,i11,lt,by)
    i7,i10 = min_max(i7,i10,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)
    i12,i14 = min_max(i12,i14,lt,by)
    i13,i15 = min_max(i13,i15,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    i10,i11 = min_max(i10,i11,lt,by)
    i12,i13 = min_max(i12,i13,lt,by)
    i14,i15 = min_max(i14,i15,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i8 = min_max(i5,i8,lt,by)
    i9,i12 = min_max(i9,i12,lt,by)
    i13,i14 = min_max(i13,i14,lt,by)

    i5,i6 = min_max(i5,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)
    i9,i10 = min_max(i9,i10,lt,by)
    i11,i12 = min_max(i11,i12,lt,by)

    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)
    i10,i11 = min_max(i10,i11,lt,by)
    i12,i13 = min_max(i12,i13,lt,by)

    i7,i8 = min_max(i7,i8,lt,by)
    i9,i10 = min_max(i9,i10,lt,by)
    return i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16
end

swapsort(x::NTuple{16,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16],lt=lt,by=by)


#Compares: 1 Depth 1 
function swapsort(i1, i2;lt=isless,by=identity)
    i1,i2 = min_max(i1,i2,lt,by)
    return i1, i2
end

swapsort(x::NTuple{2,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2],lt=lt,by=by)


#Compares: 3 Depth 3 
function swapsort(i1, i2, i3;lt=isless,by=identity)
    i1,i3 = min_max(i1,i3,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    return i1, i2, i3
end

swapsort(x::NTuple{3,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3],lt=lt,by=by)


#Compares: 5 Depth 3 
function swapsort(i1, i2, i3, i4;lt=isless,by=identity)
    i1,i4 = min_max(i1,i4,lt,by)
    i2,i3 = min_max(i2,i3,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i4 = min_max(i3,i4,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    return i1, i2, i3, i4
end

swapsort(x::NTuple{4,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4],lt=lt,by=by)


#Compares: 9 Depth 5 
function swapsort(i1, i2, i3, i4, i5;lt=isless,by=identity)
    i1,i5 = min_max(i1,i5,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i2,i5 = min_max(i2,i5,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i4 = min_max(i3,i4,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    return i1, i2, i3, i4, i5
end

swapsort(x::NTuple{5,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5],lt=lt,by=by)


#Compares: 12 Depth 5 
function swapsort(i1, i2, i3, i4, i5, i6;lt=isless,by=identity)
    i1,i6 = min_max(i1,i6,lt,by)
    i2,i4 = min_max(i2,i4,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i2,i5 = min_max(i2,i5,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i4 = min_max(i3,i4,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    return i1, i2, i3, i4, i5, i6
end

swapsort(x::NTuple{6,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6],lt=lt,by=by)


#Compares: 16 Depth 6 
function swapsort(i1, i2, i3, i4, i5, i6, i7;lt=isless,by=identity)
    i1,i7 = min_max(i1,i7,lt,by)
    i2,i6 = min_max(i2,i6,lt,by)
    i3,i4 = min_max(i3,i4,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i2,i5 = min_max(i2,i5,lt,by)
    i4,i7 = min_max(i4,i7,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    return i1, i2, i3, i4, i5, i6, i7
end

swapsort(x::NTuple{7,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7],lt=lt,by=by)


#Compares: 19 Depth 6 
function swapsort(i1, i2, i3, i4, i5, i6, i7, i8;lt=isless,by=identity)
    i1,i8 = min_max(i1,i8,lt,by)
    i2,i7 = min_max(i2,i7,lt,by)
    i3,i6 = min_max(i3,i6,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i2,i4 = min_max(i2,i4,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)
    i6,i8 = min_max(i6,i8,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    return i1, i2, i3, i4, i5, i6, i7, i8
end

swapsort(x::NTuple{8,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8],lt=lt,by=by)


#Compares: 25 Depth 7 
function swapsort(i1, i2, i3, i4, i5, i6, i7, i8, i9;lt=isless,by=identity)
    i1,i9 = min_max(i1,i9,lt,by)
    i2,i7 = min_max(i2,i7,lt,by)
    i3,i6 = min_max(i3,i6,lt,by)
    i5,i8 = min_max(i5,i8,lt,by)

    i1,i5 = min_max(i1,i5,lt,by)
    i3,i7 = min_max(i3,i7,lt,by)
    i4,i8 = min_max(i4,i8,lt,by)
    i6,i9 = min_max(i6,i9,lt,by)

    i1,i3 = min_max(i1,i3,lt,by)
    i2,i6 = min_max(i2,i6,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    i7,i9 = min_max(i7,i9,lt,by)

    i2,i4 = min_max(i2,i4,lt,by)
    i5,i7 = min_max(i5,i7,lt,by)
    i6,i8 = min_max(i6,i8,lt,by)

    i1,i2 = min_max(i1,i2,lt,by)
    i3,i5 = min_max(i3,i5,lt,by)
    i4,i6 = min_max(i4,i6,lt,by)
    i8,i9 = min_max(i8,i9,lt,by)

    i3,i4 = min_max(i3,i4,lt,by)
    i5,i6 = min_max(i5,i6,lt,by)
    i7,i8 = min_max(i7,i8,lt,by)

    i2,i3 = min_max(i2,i3,lt,by)
    i4,i5 = min_max(i4,i5,lt,by)
    i6,i7 = min_max(i6,i7,lt,by)
    return i1, i2, i3, i4, i5, i6, i7, i8, i9
end

swapsort(x::NTuple{9,T};lt=isless,by=identity) where {T} = swapsort(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9],lt=lt,by=by)


