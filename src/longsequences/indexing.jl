###
### LongSequence specific specializations of src/biosequence/indexing.jl
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# assumes `i` is positive and `bitsof(A)` is a power of 2

function Base.getindex(seq::LongSequence, part::UnitRange{<:Integer})
	@boundscheck checkbounds(seq, part)
	newseq = typeof(seq)(length(part))
	return copyto!(newseq, 1, seq, first(part), length(part))
end

function Base.setindex!(seq::LongSequence{A},
                        other::LongSequence{A},
                        locs::AbstractVector{<:Integer}) where {A}
    @boundscheck checkbounds(seq, locs)
    checkdimension(other, locs)
    return unsafe_setindex!(seq, other, locs)
end

function Base.setindex!(seq::LongSequence{A},
                        other::LongSequence{A},
                        locs::AbstractVector{Bool}) where {A}
    @boundscheck checkbounds(seq, locs)
    checkdimension(other, locs)
    return unsafe_setindex!(seq, other, locs)
end

function Base.setindex!(seq::LongSequence{A},
                        other::LongSequence{A},
                        locs::UnitRange{<:Integer}) where {A}
    @boundscheck checkbounds(seq, locs)
    checkdimension(other, locs)
    return copyto!(seq, locs.start, other, 1, length(locs))
end

function Base.setindex!(seq::LongSequence{A},
			            other::LongSequence{A}, ::Colon) where {A}
    return setindex!(seq, other, 1:lastindex(seq))
end

function Base.setindex!(seq::LongSequence{A}, x, ::Colon) where {A}
    return setindex!(seq, x, 1:lastindex(seq))
end

@inline function encoded_setindex!(s::LongSequence, v::Unsigned, i::BitIndex)
    vi, off = i
    data = encoded_data(s)
    bits = @inbounds data[vi]
	v_ = v % encoded_data_eltype(s)
    @inbounds data[vi] = (v_ << off) | (bits & ~(bindata_mask(s) << off))
    return s
end

@inline function unsafe_setindex!(seq::LongSequence, x, locs::AbstractVector{<:Integer})
	bin = encode(Alphabet(seq), convert(eltype(seq), x))
    for i in locs
        encoded_setindex!(seq, bin, bitindex(seq, i))
    end
    return seq
end

function unsafe_setindex!(seq::LongSequence, x, locs::AbstractVector{Bool})
    bin = encode(Alphabet(seq), convert(eltype(seq), x))
    i = j = 0
    while true
        i = findnext(locs, i + 1)
        if i === nothing
            break
        end
        encoded_setindex!(seq, bin, i)
    end
    return seq
end

function unsafe_setindex!(seq::LongSequence{A}, other::LongSequence{A},
	                      locs::AbstractVector{<:Integer}) where {A <: Alphabet}
    for (i, x) in zip(locs, other)
        unsafe_setindex!(seq, x, i)
    end
    return seq
end

function unsafe_setindex!(seq::LongSequence{A},
                          other::LongSequence{A},
                          locs::AbstractVector{Bool}) where {A <: Alphabet}
    i = j = 0
    while true
        i = findnext(locs, i + 1)
        if i === nothing
            break
        end
        unsafe_setindex!(seq, other[j += 1], i)
    end
    return seq
end

@inline function encoded_setindex!(seq::LongSequence{A},
				   bin::Unsigned, i::Integer) where {A <: Alphabet}
    j, r = bitindex(seq, i)
	bin_ = bin % encoded_data_eltype(seq)
    data = encoded_data(seq)
    @inbounds data[j] = (bin_ << r) | (data[j] & ~(bindata_mask(seq) << r))
    return seq
end
