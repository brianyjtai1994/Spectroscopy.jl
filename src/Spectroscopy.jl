module Spectroscopy

const VecI  = AbstractVector # Input  Vector
const VecO  = AbstractVector # Output Vector
const VecB  = AbstractVector # Buffer Vector
const VecIO = AbstractVector # In/Out Vector
const MatI  = AbstractMatrix # Input  Matrix
const MatO  = AbstractMatrix # Output Matrix
const MatB  = AbstractMatrix # Buffer Matrix
const MatIO = AbstractMatrix # In/Out Matrix

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

include("./TimeDelay.jl")
include("./DifferentialData.jl")

end # module
