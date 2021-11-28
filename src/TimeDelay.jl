export TimeDelay, ps2μm!, μm2ps!

ps2μm(t::Real) = round(149.896225 * t)
μm2ps(d::Real) = d / 149.896225

for f in (:ps2μm, :μm2ps)
    @eval begin
        function $(Symbol(f, !))(arr::VecI)
            @simd for i in eachindex(arr)
                @inbounds arr[i] = $f(arr[i])
            end
            return nothing
        end
    end
end

struct TimeDelay <: TimeAxis
    dat::Vector{Float64}

    function TimeDelay(src::VecI, unit::String="ps")
        dat = similar(src)
        if unit ≡ "ps"
            @simd for i in eachindex(dat)
                @inbounds dat[i] = ps2μm(src[i])
            end
        elseif unit ≡ "μm"
            @simd for i in eachindex(dat)
                @inbounds dat[i] = src[i]
            end
        else
            error("TimeDelay(..., unit = $unit) is not allowed.")
        end
        return new(dat)
    end
end
