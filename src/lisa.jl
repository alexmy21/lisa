# This module is greatly inspired by the implementation of the HLL algorithm in the Julia library:
# https://github.com/jakobnissen/Probably.jl/blob/master/src/hyperloglog/hyperloglog.jl

# Below is the header from this file
# From Flajolet, Philippe; Fusy, Éric; Gandouet, Olivier; Meunier, Frédéric (2007)
# DOI: 10.1.1.76.4286
# With algorithm improvements by Google (https://ai.google/research/pubs/pub40671)

# Principle:
# When observing N distinct uniformly distributed integers, the expected maximal
# number of leading zeros in the integers is log(2, N), with large variation.
# To cut variation, we keep 2^P counters, each keeping track of N/2^P
# observations. The estimated N for each counter is averaged using harmonic mean.
# Last, corrections for systematic bias are added, one multiplicative and one
# additive factor.
# To make the observations uniformly distributed integers, we hash them.

# We made sugnificant changes to the original implementation:
# - We use a BitVector instead of a UInt8 for the counters
# - We implemented additional operators to support set operations, like union (union), intersection(intersect), difference(diff), 
#   and equality (isequal). Now they work the same way as they work for sets
# - We also renamed some of the original operators to be more consistent with HyperLogLog terminology


include("constants.jl")

    struct HllSet{P}
        counts::Vector{BitVector}

        function HllSet{P}() where {P}
            validate_P(P)
            n = calculate_n(P)
            counts = create_bitsets(n)
            return new(counts)
        end

        function validate_P(P)
            isa(P, Integer) || throw(ArgumentError("P must be integer"))
            (P < 4 || P > 18) && throw(ArgumentError("P must be between 4 and 18"))
        end

        function calculate_n(P)
            return 1 << P
        end

        function create_bitsets(n)
            return [falses(64) for _ in 1:n]
        end
    end

    HllSet{P}(::Type{T}) where {P, T} = HllSet{10}()

    Base.show(io::IO, x::HllSet{P}) where {P} = print(io, "HllSet{$(P)}()")

    Base.sizeof(::Type{HllSet{P}}) where {P} = 1 << P
    Base.sizeof(x::HllSet{P}) where {P} = sizeof(typeof(x))

    function Base.union!(dest::HllSet{P}, src::HllSet{P}) where {P}
        length(dest.counts) == length(src.counts) || throw(ArgumentError("HllSet{P} must have same size"))
        for i in 1:length(dest.counts)
            dest.counts[i] = dest.counts[i] .| src.counts[i]
        end
        return dest
    end

    function Base.copy!(dest::HllSet{P}, src::HllSet{P}) where {P}
        length(dest.counts) == length(src.counts) || throw(ArgumentError("HllSet{P} must have same size"))
        for i in 1:length(dest.counts)
            dest.counts[i] = src.counts[i]
        end
        return dest
    end

    function Base.union(x::HllSet{P}, y::HllSet{P}) where {P} 
        length(x.counts) == length(y.counts) || throw(ArgumentError("HllSet{P} must have same size"))
        z = HllSet{P}()
        for i in 1:length(x.counts)
            z.counts[i] = x.counts[i] .| y.counts[i]
        end
        return z
    end

    function Base.intersect(x::HllSet{P}, y::HllSet{P}) where {P} 
        length(x.counts) == length(y.counts) || throw(ArgumentError("HllSet{P} must have same size"))
        z = HllSet{P}()
        for i in 1:length(x.counts)
            z.counts[i] = x.counts[i] .& y.counts[i]
        end
        return z
    end

    function Base.diff(x::HllSet{P}, y::HllSet{P}) where {P} 
        length(x.counts) == length(y.counts) || throw(ArgumentError("HllSet{P} must have same size"))
        z = HllSet{P}()
        for i in 1:length(x.counts)
            z.counts[i] = x.counts[i] .& .~(y.counts[i])
        end
        return z
    end

    function Base.isequal(x::HllSet{P}, y::HllSet{P}) where {P} 
        length(x.counts) == length(y.counts) || throw(ArgumentError("HllSet{P} must have same size"))
        for i in 1:length(x.counts)
            x.counts[i] == y.counts[i] || return false
        end
        return true
    end

    Base.isempty(x::HllSet{P}) where {P} = all(x -> x == 0,  x.counts)

    function getbin(hll::HllSet{P}, x::UInt) where {P} 
        x = x >>> (8 * sizeof(UInt) - P) + 1
        str = replace(string(x, base = 16), "0x" => "")
        return parse(Int, str, base = 16)
    end

    function getzeros(hll::HllSet{P}, x::UInt) where {P}
        or_mask = ((UInt(1) << P) - 1) << (8 * sizeof(UInt) - P)
        return trailing_zeros(x | or_mask) + 1
    end

    function Base.push!(hll::HllSet{P}, x::Any) where {P}
        h = hash(x)
        bin = getbin(hll, h)
        idx = getzeros(hll, h)
        hll.counts[bin][idx] = true
        return hll
    end

    function Base.push!(hll::HllSet{P}, values...) where {P}
        for value in values
            push!(hll, value)
        end
        return hll
    end
        
    α(x::HllSet{P}) where {P} =
        if P == 4
            return 0.673
        elseif P == 5
            return 0.697
        elseif P == 6
            return 0.709
        else
            return 0.7213 / (1 + 1.079 / sizeof(x))
        end

    
    function bias(::HllSet{P}, biased_estimate) where {P}
        # For safety - this is also enforced in the HLL constructor
        if P < 4 || P > 18
            error("We only have bias estimates for P ∈ 4:18")
        end
        rawarray = @inbounds RAW_ARRAYS[P - 3]
        biasarray = @inbounds BIAS_ARRAYS[P - 3]
        firstindex = searchsortedfirst(rawarray, biased_estimate)
        # Raw count large, no need for bias correction
        if firstindex == length(rawarray) + 1
            return 0.0
            # Raw count too small, cannot be corrected. Maybe raise error?
        elseif firstindex == 1
            return @inbounds biasarray[1]
            # Else linearly approximate the right value for bias
        else
            x1, x2 = @inbounds rawarray[firstindex - 1], @inbounds rawarray[firstindex]
            y1, y2 = @inbounds biasarray[firstindex - 1], @inbounds biasarray[firstindex]
            delta = @fastmath (biased_estimate - x1) / (x2 - x1) # relative distance of raw from x1
            return y1 + delta * (y2 - y1)
        end
    end

    function maxidx(vec::BitVector)
        
        for i in length(vec):-1:1
            if vec[i]
                return i
            end
        end
        return 0
    end

    function Base.count(x::HllSet{P}) where {P}
        # Harmonic mean estimates cardinality per bin. There are 2^P bins
        harmonic_mean = sizeof(x) / sum(1 / 1 << maxidx(i) for i in x.counts)
        biased_estimate = α(x) * sizeof(x) * harmonic_mean
        return round(Int, biased_estimate - bias(x, biased_estimate))
    end

    greet() = print("Hello World!")