###
### Mer specific specializations of src/biosequence/predicates.jl
###

Base.cmp(x::T, y::T) where {T<:Kmer} = cmp(packed_data(x), packed_data(y))
Base.:(==)(x::T, y::T) where {T<:Kmer} = packed_data(x) == packed_data(y)
Base.isless(x::T, y::T) where {T<:Kmer} = isless(packed_data(x), packed_data(y))

function Base.hash(x::Kmer{A,K,N}, h::UInt) where {A<:NucleicAcidAlphabet{2},K,N}
    return Base.hash(packed_data(x) ⊻ K ⊻ h)
end
