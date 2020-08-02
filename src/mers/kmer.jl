###
### Type definition
###

"""
    Kmer{A<:NucleicAcidAlphabet{2},K,N} <: BioSequence{A}
A parametric, immutable, bitstype for representing Kmers - short sequences.
Given the number of Kmers generated from raw sequencing reads, avoiding
repetetive memory allocation and triggering of garbage collection is important,
as is the ability to effectively pack Kmers into arrays and similar collections.
In julia this means an immutable bitstype must represent such shorter Kmer
sequences. Thankfully this is not much of a limitation - kmers are rarely
manipulated and so by and large don't have to be mutable like `LongSequence`s.
Excepting their immutability, they fulfill the rest of the API and behaviours
expected from a concrete `BioSequence` type, and non-mutating transformations
of the type are still defined.
!!! warning
    Given their immutability, `setindex` and mutating sequence transformations
    are not implemented for kmers e.g. `reverse_complement!`. 
!!! tip
    Note that some sequence transformations that are not mutating are
    available, since they can return a new kmer value as a result e.g.
    `reverse_complement`. 
"""
struct Kmer{A<:NucleicAcidAlphabet{2},K,N} <: BioSequence{A}
    data::NTuple{N,UInt64}
end


# Aliases
"Shortcut for the type `Kmer{DNAAlphabet{2},K,N}`"
const DNAKmer{K,N} = Kmer{DNAAlphabet{2},K,N}

"Shortcut for the type `DNAKmer{27,1}`"
const DNAKmer27 = DNAKmer{27,1}

"Shortcut for the type `DNAKmer{31,1}`"
const DNAKmer31 = DNAKmer{31,1}

"Shortcut for the type `DNAKmer{63,2}`"
const DNAKmer63 = DNAKmer{63,2}



###
### Base Functions
###

"""
    kmertype(::Type{Kmer{A,K}}) where {A,K}
Resolve and incomplete kmer typing, computing the N parameter of
`Kmer{A,K,N}`, given only `Kmer{A,K}`.
## Example
```julia
julia> DNAKmer{63}
Kmer{DNAAlphabet{2},63,N} where N
julia> kmertype(DNAKmer{63})
Kmer{DNAAlphabet{2},63,2}
```
"""
@inline function kmertype(::Type{Kmer{A,K}}) where {A,K}
    return Kmer{A,K,ifelse(rem(2K, 64) != 0, div(2K, 64) + 1, div(2K, 64))}
end

@inline capacity(::Type{Kmer{A,K,N}}) where {A,K,N} = div(64N, 2)
@inline capacity(seq::Kmer) = capacity(typeof(seq))
@inline n_unused(::Type{Kmer{A,K,N}}) where {A,K,N} = capacity(Kmer{A,K,N}) - K
@inline n_unused(seq::Kmer) = n_unused(typeof(seq))

@inline function checkmer(::Type{Kmer{A,K,N}}) where {A,K,N}
    n = ifelse(rem(2K, 64) != 0, div(2K, 64) + 1, div(2K, 64))
    if n !== N
        # This has been significantly changed conceptually from before. Now we
        # don't just check K, but enforce the most appropriate N for K.
        throw(ArgumentError("Bad kmer parameterisation. For K = $K, N should be $n"))
    end
end

#@inline Base.eltype(::Type{Kmer{A,K,N}}) where {A,K,N} = eltype(A)
@inline Base.length(x::Kmer{A,K,N}) where {A,K,N} = K
@inline Base.summary(x::Kmer{DNAAlphabet{2},K,N}) where {K,N} = string("DNA ", K, "-mer")
@inline Base.summary(x::Kmer{RNAAlphabet{2},K,N}) where {K,N} = string("RNA ", K, "-mer")

function Base.typemin(::Type{Kmer{A,K,N}}) where {A,K,N}
    checkmer(Kmer{A,K,N})
    return Kmer{A,K,N}(ntuple(i -> zero(UInt64), N))
end

function Base.typemax(::Type{Kmer{A,K,N}}) where {A,K,N}
    checkmer(Kmer{A,K,N})
    return Kmer{A,K,N}((typemax(UInt64) >> (64N - 2K), ntuple(i -> typemax(UInt64), N - 1)...))
end

function Base.rand(::Type{Kmer{A,K,N}}) where {A,K,N}
    checkmer(Kmer{A,K,N})
    return Kmer{A,K,N}((rand(UInt64) >> (64N - 2K), ntuple(i -> rand(UInt64), N - 1)...))
end

function Base.rand(::Type{T}, size::Integer) where {T<:Kmer}
    return [rand(T) for _ in 1:size]
end

###
### Old Mer Base Functions - not transferred to new type.
###
#@inline encoded_data_type(::Type{Mer{A,K}}) where {A,K} = UInt64
#@inline encoded_data_type(::Type{BigMer{A,K}}) where {A,K} = UInt128
#@inline encoded_data_type(x::AbstractMer) = encoded_data_type(typeof(x))
#@inline encoded_data(x::AbstractMer) = reinterpret(encoded_data_type(typeof(x)), x)
#@inline ksize(::Type{T}) where {A,K,T<:AbstractMer{A,K}} = K
#@inline Base.unsigned(x::AbstractMer) = encoded_data(x)
#Base.:-(x::AbstractMer, y::Integer) = typeof(x)(encoded_data(x) - y % encoded_data_type(x))
#Base.:+(x::AbstractMer, y::Integer) = typeof(x)(encoded_data(x) + y % encoded_data_type(x))
#Base.:+(x::AbstractMer, y::AbstractMer) = y + x
#Alphabet(::Type{Mer{A,K} where A<:NucleicAcidAlphabet{2}}) where {K} = Any

###
### Constructors
###

# Create a Mer from a sequence.
function (::Type{Kmer{A,K,N}})(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet{2},K,N}
    seqlen = length(seq)
    if seqlen != K
        throw(ArgumentError("seq does not contain the correct number of nucleotides ($seqlen ≠ $K)"))
    end
    checkmer(Kmer{A,K,N})

    # Construct the head.
    bases_in_head = div(64 - (64N - 2K), 2)
    head = zero(UInt64)
    @inbounds for i in 1:bases_in_head
        nt = convert(eltype(typeof(seq)), seq[i])
        head = (head << 2) | UInt64(twobitnucs[reinterpret(UInt8, nt) + 0x01])
    end
    
    # And the rest of the sequence
    idx = Ref(bases_in_head + 1)
    
    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for i in 1:32
            nt = convert(eltype(typeof(seq)), seq[idx[]])
            body = (body << 2) | UInt64(twobitnucs[reinterpret(UInt8, nt) + 0x01])
            idx[] += 1
        end
        return body
    end

    return Kmer{A,K,N}((head, tail...))
end

@inline function inbounds_getindex(x::Kmer{A,K,N}, i::Integer) where {A,K,N}
    # Emulation of BitIndex type
    i′ = i + div(64N - 2K, 2)
    val = (i′ - 1) << 1
    idx = (val >> 6) + 1
    off = 62 - (val & (UInt8(64) - 0x01))
    @inbounds begin
        chunk = x.data[idx]
    end
    bits = (chunk >> off) & UInt64(3)
    return reinterpret(eltype(x), 0x01 << bits)
end

include("predicates.jl")
#include("counting.jl") # TODO: Update for NTuple based kmers.
include("transformations.jl")


###
### Kmer de-bruijn neighbors
###

# TODO: Decide on this vs. old iterator pattern. I like the terseness of the code vs defining an iterator. Neither allocate.
fw_neighbors(kmer::Kmer) = ntuple(Val{4}(), i -> kmer << ACGT[i])


#=
# Neighbors on a de Bruijn graph
struct KmerNeighborIterator{S<:Kmer}
    x::S
end

"""
    neighbors(kmer::S) where {S<:Kmer}

Return an iterator through skip-mers neighboring `skipmer` on a de Bruijn graph.
"""
neighbors(kmer::Kmer) = KmerNeighborIterator{typeof(kmer)}(kmer)

Base.length(::KmerNeighborIterator) = 4
Base.eltype(::Type{KmerNeighborIterator{S}}) where {S<:Kmer} = S

function Base.iterate(it::KmerNeighborIterator{S}, i::UInt64 = 0) where {S<:Kmer}
    if i == 4
        return nothing
    else
        #return S((encoded_data(it.x) << 2) | i), i + 1
        return it.x << 1, i + one(UInt64)
    end
end
=#

###
### String literals
###

macro mer_str(seq, flag)
    seq′ = remove_newlines(seq)
    if flag == "dna" || flag == "d"
        T = kmertype(DNAKmer{length(seq′)})
        return T(seq′)
    elseif flag == "rna" || flag == "r"
        T = kmertype(RNAKmer{length(seq′)})
        return T(seq′)
    else
        error("Invalid type flag: '$(flag)'")
    end
end

macro mer_str(seq)
    seq′ = remove_newlines(seq)
    T = kmertype(DNAKmer{length(seq′)})
    return T(seq′)
end
