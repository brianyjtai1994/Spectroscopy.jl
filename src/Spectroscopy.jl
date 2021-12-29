module Spectroscopy

const VecI  = AbstractVector # Input  Vector
const VecO  = AbstractVector # Output Vector
const VecB  = AbstractVector # Buffer Vector
const VecIO = AbstractVector # In/Out Vector
const MatI  = AbstractMatrix # Input  Matrix
const MatO  = AbstractMatrix # Output Matrix
const MatB  = AbstractMatrix # Buffer Matrix
const MatIO = AbstractMatrix # In/Out Matrix

using DelimitedFiles

abstract type TimeResolvedData             end
abstract type TimeAxis <: TimeResolvedData end

macro def(genre::Symbol, ex::Union{Expr,Symbol}, vars::Symbol...)
    o = genre == :prop ? :(::) : genre == :vars ? :(=) : error("@def(genre = $genre, ...) is invalid.")
    n = length(vars)
    e = Vector{Expr}(undef, n)
    @inbounds for i in 1:n
        e[i] = Expr(o, vars[i], ex)
    end
    return Expr(:escape, Expr(:block, e...))
end

macro get(obj::Symbol, vars::Symbol...)
    n = length(vars)
    e = Vector{Expr}(undef, n)
    @inbounds for i in 1:n
        vari = vars[i]
        e[i] = :($vari = $obj.$vari)
    end
    return Expr(:escape, Expr(:block, e...))
end

function polyinterp(x::Real, xv::VecI, yv::VecI, bv::VecB, n::Int)
    one2n = eachindex(1:n)
    @inbounds for i in one2n
        x == xv[i] && return yv[i]
    end
    Δx = Inf
    yp = 0.0
    @inbounds for i in one2n
        δx = abs(x - xv[i])
        δx < Δx && (Δx = δx; yp = yv[i])
    end
    @simd for i in one2n
        @inbounds bv[i] = yv[i] - yp
    end
    @inbounds for k in 1:n-1, i in 1:n-k
        bv[i] += (bv[i] - bv[i+1]) * (x - xv[i]) / (xv[i] - xv[i+k])
    end
    @inbounds return yp + bv[1]
end

function polyinterp!(ydes::VecO{Ty}, xdes::VecO{Tx}, ysrc::VecI, xsrc::VecI; order::Int=3) where {Tx<:Real,Ty<:Real}
    lx = 1                     # left index
    rx = sz = order + 1        # right index
    ia = sz >> 1 + 1           # anchor index
    mx = size(xsrc, 1)         # upper bound of window
    xb = Vector{Tx}(undef, sz) # xsrc buffer
    yb = Vector{Ty}(undef, sz) # ysrc buffer
    tb = Vector{Ty}(undef, sz) # interpolation buffer
    @inbounds for i in 1:sz
        xb[i] = xsrc[i]
        yb[i] = ysrc[i]
    end
    @inbounds xanchor = xb[ia]
    @inbounds for i in eachindex(xdes)
        xi = xdes[i]
        #### moving interpolation window
        while xi > xanchor && rx < mx
            lx += 1
            rx += 1
            ix  = 1
            jx  = 2
            while ix < sz
                @inbounds xb[ix] = xb[jx]
                @inbounds yb[ix] = yb[jx]
                ix += 1
                jx += 1
            end
            @inbounds xb[ix]  = xsrc[rx]
            @inbounds yb[ix]  = ysrc[rx]
            @inbounds xanchor = xb[ia]
        end
        ydes[i] = polyinterp(xi, xb, yb, tb, sz)
    end
    return nothing
end

export polyinterp!

include("./TimeDelay.jl")
include("./DifferentialData.jl")

end # module
