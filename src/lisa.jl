include("constants.jl")

module lisa

    struct HllSet{P}
        counts::Vector{BitSet}

        function HllSet{P}() where {P}
            isa(P, Integer) || throw(ArgumentError("P must be integer"))
            (P < 4 || P > 18) && throw(ArgumentError("P must be between 4 and 18"))
            return new(zeros(UInt8, sizeof(new([BitSet() for _ in 1:HllSet{P}]))))
        end
    end

    HllSet{P}(::Type{T}) where {P, T} = HllSet{10}()

    Base.show(io::IO, hll::HllSet) = print(io, "HllSet{", P, "}(", length(hll.counts), ")")
    Base.sizeof(::Type{HllSet{P}}) where {P}= 1 << P
    Base.sizeof(hll::HllSet) = sizeof(typeof(hll))

    function Base.update!(dest::HllSet{P}, src::HllSet{P}) where {P}
        length(dest.counts) == length(src.counts) || throw(ArgumentError("HllSet{P} must have same size"))
        for i in 1:length(dest.counts)
            dest.counts[i] |= src.counts[i]
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
            z.counts[i] = x.counts[i] | y.counts[i]
        end
        return z
    end

    function Base.intersection(x::HllSet{P}, y::HllSet{P}) where {P} 
        length(x.counts) == length(y.counts) || throw(ArgumentError("HllSet{P} must have same size"))
        z = HllSet{P}()
        for i in 1:length(x.counts)
            z.counts[i] = x.counts[i] & y.counts[i]
        end
        return z
    end

    function Base.setdiff(x::HllSet{P}, y::HllSet{P}) where {P} 
        length(x.counts) == length(y.counts) || throw(ArgumentError("HllSet{P} must have same size"))
        z = HllSet{P}()
        for i in 1:length(x.counts)
            z.counts[i] = x.counts[i] & ~y.counts[i]
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

    function Base.hash(x::Any)  
        hash(x)
    end

    getbin(hll::HllSet{P}, x::UInt) where {P} = x >>> (8 * sizeof(UInt) - P) + 1

    function getzeros(hll::HllSet{P}, x::UInt) where {P}
        or_mask = ((UInt(1) << P) - 1) << (8 * sizeof(UInt) - P)
        return trailing_zeros(x | or_mask) + 1
    end

    function Base.push!(hll::HllSet, x)
        h = hash(x)
        bin = getbin(hll, h)
        @inbounds hll.counts[bin] |= BitSet(1 << getzeros(hll, h))
        return hll
    end

    function Base.push!(hll::HllSet, values...)
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

    function Base.length(x::HllSet{P}) where {P}
        # Harmonic mean estimates cardinality per bin. There are 2^P bins
        harmonic_mean = sizeof(x) / sum(1 / 1 << maximum(i) for i in x.counts)
        biased_estimate = α(x) * sizeof(x) * harmonic_mean
        return round(Int, biased_estimate - bias(x, biased_estimate))
    end

    greet() = print("Hello World!")

end # module lisa