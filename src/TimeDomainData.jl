export TimeDomainData, ComplexSpectrum

struct TimeDomainData <: TimeResolvedData
    td::Vector{Float64} # time delay
    ΔE::Vector{Float64} # electric field difference
    Δt::Float64

    function TimeDomainData(fname::String; args...)
        dat = readdlm(fname; args...)
        td  = collect(range(dat[1,1], dat[end,1], length=clp2(size(dat, 1))))
        ΔE  = similar(td)
        polyinterp!(ΔE, td, view(dat, :, 2), view(dat, :, 1))
        return new(td, ΔE, td[2] - td[1])
    end
end

struct ComplexSpectrum <: TimeResolvedData
    @def prop Vector{Float64} fs Re Im Rs θs; sz::Int

    function ComplexSpectrum(ref::NamedTuple{(:frq, :cxy), Tuple{Vector{Float64}, Matrix{Float64}}},
                             sam::NamedTuple{(:frq, :cxy), Tuple{Vector{Float64}, Matrix{Float64}}};
                             fmin::Real=0.1, fmax::Real=3.0)
        frq = ref.frq
        length(frq) ≠ length(sam.frq) && error("length(ref.frq) ≠ length(sam.frq)")

        lx = length(frq) >> 1
        @inbounds while frq[lx] < fmin
            lx += 1
        end
        rx = lx
        @inbounds while !(frq[rx] > fmax)
            rx += 1
        end
        sz = rx - lx

        @def vars Vector{Float64}(undef, sz) fs Re Im Rs θs

        @simd for i in eachindex(fs)
            @inbounds fs[i] = frq[lx+i]
        end

        rxy = ref.cxy
        sxy = sam.cxy

        @inbounds for i in eachindex(fs)
            Re[i], Im[i] = cdiv(sxy[lx+i,1], sxy[lx+i,2], rxy[lx+i,1], rxy[lx+i,2])
        end

        #### Phases unwrapping
        halfπ = 0.5 * π
        @inbounds θs[1] = prevθ = atan(Im[1], Re[1])
        @inbounds Rs[1] = 0.0
        @inbounds for i in 2:sz
            θs[i] = thisθ = atan(Im[i], Re[i])
            diff = thisθ - prevθ
            θmod = rem2pi(diff + halfπ, RoundDown) - halfπ
            diff > 0.0 && θmod == -halfπ && (θmod = halfπ)
            Rs[i] = abs(diff) < halfπ ? Rs[i-1] : θmod - diff + Rs[i-1]
            prevθ = thisθ
        end

        @simd for i in eachindex(fs)
            @inbounds θs[i] += Rs[i]
        end

        @simd for i in eachindex(fs)
            @inbounds Rs[i] = apy2(Re[i], Im[i])
        end

        return new(fs, Re, Im, Rs, θs, sz)
    end
end
