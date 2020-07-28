###
### Internal NTuple manipulation methods
###

"""
    _cliphead(by::Integer, head::UInt64, tail...)

A method used to mask the first `by` MSB's in `head`, before catting it with
tail to return a NTuple.

This is used internally to mask the first `by` bits in the first word of a 
NTuple of UInt64's.

Notably it's used when constructing a Kmer from an existing NTuple of UInt64
"""
@inline function _cliphead(by::Integer, head::UInt64, tail...)
    return (head & (typemax(UInt64) >> by), tail...)
end

@inline function choptail(x::NTuple{N,UInt64}) where {N}
    ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        return @inbounds x[i]
    end
end

@inline function setlast(x::NTuple{N,UInt64}, nt::DNA) where {N}
    @inbounds begin
        bits = UInt64(twobitnucs[reinterpret(UInt8, nt) + 0x01])
        tail = (x[N] & (typemax(UInt64) - UInt64(3))) | bits
    end
    return (choptail(x)..., tail)
end

# Bit-parallel element nucleotide complementation
@inline function complement_bitpar(x::NTuple{N,U}, a::A) where {N,U<:Unsigned,A<:NucleicAcidAlphabet{2}}
    return _complement_bitpar(a, x...)
end

@inline function _complement_bitpar(a::A, head::UInt64, tail...) where {A<:NucleicAcidAlphabet{2}}
    
    return (~head, _complement_bitpar(a, tail...)...)
end

@inline _complement_bitpar(a::A) where {A<:NucleicAcidAlphabet{2}} = ()

@inline function reverse(bpe::BitsPerElem{B}, x::NTuple{N,UInt64}) where {B,N}
    return _reverse(bpe, x...)
end

@inline function _reverse(bpe::BitsPerElem{N}, head::UInt64, tail...) where {N}
    return (_reverse(bpe, tail...)..., reversebits(head, bpe))
end

@inline _reverse(::BitsPerElem{N}) = ()

#=
@inline function reversebits(x::T, ::BitsPerElem{2}) where {T<:Base.BitUnsigned}
     mask = 0x33333333333333333333333333333333 % T
     x = ((x >> 2) & mask) | ((x & mask) << 2)
     return reversebits(x, BitsPerElem{4}())
end

@inline function reversebits(x::T, ::BitsPerElem{4}) where {T<:Base.BitUnsigned}
     mask = 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F % T
     x = ((x >> 4) & mask) | ((x & mask) << 4)
     return bswap(x)
end
=#

@inline function rightshift_carry(x::NTuple{N,UInt64}, nbits::Integer) where {N}
    return _rightshift(n, zero(UInt64), x...)
end

@inline function _rightshift_carry(nbits::Integer, carry::UInt64, head::UInt64, tail...)
    return ((head >> nbits) | carry, _rightshift(nbits, (head & ((one(UInt64) << nbits) - 1)) << (64 - nbits), tail...)...)
end

@inline _rightshift_carry(nbits::Integer, carry::UInt64) = ()


# TODO: Work this into a left version of rightshift_carry... which it already is, but only shifts 2 bits, no `nbits` parameter is currently accepted.
"""
    shiftleft

It is important to be able to efficiently shift the all the nucleotides in a kmer
one space to the left or right, as it is a key operation in iterating through
de bruijn graph neighbours or in building kmers a nucleotide at a time.
"""
function shiftleft end

@inline function shiftleft(x::Kmer{A,K,N}) where {A,K,N}
    _, newbits = _shiftleft(x.data...)
    # TODO: The line below is a workaround for julia issues #29114 and #3608
    newbits′ = newbits isa UInt64 ? (newbits,) : newbits
    return Kmer{A,K,N}(_cliphead(64N - 2K, newbits′...))
end

@inline function shiftleft(x::NTuple{N,UInt64}) where {N}
    _, newbits = _shiftleft(x...)
    # TODO: The line below is a workaround for julia issues #29114 and #3608
    return newbits isa UInt64 ? (newbits,) : newbits
end

@inline function _shiftleft(head::UInt64, tail...)
    carry, newtail = _shiftleft(tail...)
    # TODO: The line below is a workaround for julia issues #29114 and #36087
    newtail′ = newtail isa UInt64 ? (newtail,) : newtail 
    return head >> 62, ((head << 2) | carry, newtail′...)
end

@inline _shiftleft(head::UInt64) = (head & 0xC000000000000000) >> 62, head << 2


###
### Transformation methods
###

@inline Base.:(>>)(seq::Kmer{A,K,N}, nbases::Integer) where {A,K,N} = rightshift_carry(seq.data, 2nbases)


@inline function Base.:(<<)(seq::Kmer{DNAAlphabet{2},K,N}, nuc::DNA) where {K,N}
    return Kmer{A,K,N}(_cliphead(64N - 2K, setlast(shiftleft(packed_data(seq)), nuc)...))
end

"""
    complement(seq::T) where {T<:Kmer}
Return the complement of a short sequence type `x`.
"""
@inline function complement(seq::T) where {T<:Kmer}
    return T(complement_bitpar(packed_data(seq), Alphabet(seq)))
end

"""
    reverse(seq::Kmer{A,K,N}) where {A,K,N}
Return the reverse of short sequence type variable `seq`.
"""
function Base.reverse(seq::Kmer{A,K,N}) where {A,K,N}
    bits = reverse(BitsPerElem(seq), packed_data(seq))
    return T(rightshift_carry(sbits, n_unused(seq)))
end

# TODO: Check BioSequences.jl doesen't basically give this generic def for free for any BioSequence anyway.
"""
    reverse_complement(x::Kmer)

Return the reverse complement of `x`.
"""
reverse_complement(x::Kmer) = complement(reverse(x))

"""
    canonical(x::Kmer)

Return the canonical sequence of `x`.

A canonical sequence is the numerical lesser of a k-mer and its reverse complement.
This is useful in hashing/counting sequences in data that is not strand specific,
and thus observing the short sequence is equivalent to observing its reverse complement.
"""
@inline canonical(x::Kmer) = min(x, reverse_complement(x))

###
### Old Mer specific specializations of src/biosequence/transformations.jl
### - not currently transferred to new type.

# TODO: Sort this and decide on transferring to new NTuple based kmers or no.

#=
function swap(x::T, i, j) where {T<:AbstractMer}
    i = 2 * length(x) - 2i
    j = 2 * length(x) - 2j
    b = encoded_data(x)
    x = ((b >> i) ⊻ (b >> j)) & encoded_data_type(x)(0x03)
    return T(b ⊻ ((x << i) | (x << j)))
end


function Random.shuffle(x::T) where {T<:AbstractMer}
    # Fisher-Yates shuffle for mers.
    j = lastindex(x)
    for i in firstindex(x):(j - 1)
        j′ = rand(i:j)
        x = swap(x, i, j′)
    end
    return x
end
=#