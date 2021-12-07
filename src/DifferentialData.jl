export DifferentialData, average!

struct DifferentialData <: TimeResolvedData
    @def prop Vector{Float64} real δreal imag δimag
    logs::Dict{Symbol,Int} # [dims, avg. num]

    DifferentialData(n::Int) = new(
        Vector{Float64}(undef, n), Vector{Float64}(undef, n),
        Vector{Float64}(undef, n), Vector{Float64}(undef, n),
        Dict(:dims => n, :count => 0)
    )

    function DifferentialData(s::String, c::Int=1; Ireal::Int=2)
        data = readdlm(s)
        dims = size(data, 1)
        real = Vector{Float64}(undef, dims)
        @simd for i in eachindex(1:dims)
            @inbounds real[i] = data[i, Ireal]
        end
        return new(
            real, zeros(Float64, dims), zeros(Float64, dims), zeros(Float64, dims),
            Dict(:dims => dims, :count => c)
        )
    end
end

# Perform a single step of Wolford algorithm
function welford_step(μ::Real, s::Real, v::Real, c::Real)
    isone(c) && return v, zero(v)
    s = s * (c - 1)
    m = μ + (v - μ) / c
    s = s + (v - μ) * (v - m)
    μ = m
    return μ, s / (c - 1)
end

function average!(des::DifferentialData, src::MatI, Ireal::Int, Iimag::Int, Iref::Int)
    @get des real δreal imag δimag logs

    one2n = eachindex(1:logs[:dims])
    c = logs[:count] + 1

    @inbounds for i in one2n
        ref = src[i, Iref]
        real[i], δreal[i] = welford_step(real[i], δreal[i], src[i, Ireal] / ref, c)
        imag[i], δimag[i] = welford_step(imag[i], δimag[i], src[i, Iimag] / ref, c)
    end

    logs[:count] = c
    return nothing
end
