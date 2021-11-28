module Spectroscopy

const VecI  = AbstractVector # Input  Vector
const VecO  = AbstractVector # Output Vector
const VecB  = AbstractVector # Buffer Vector
const VecIO = AbstractVector # In/Out Vector
const MatI  = AbstractMatrix # Input  Matrix
const MatO  = AbstractMatrix # Output Matrix
const MatB  = AbstractMatrix # Buffer Matrix
const MatIO = AbstractMatrix # In/Out Matrix

abstract type TimeAxis end

include("./TimeDelay.jl")

end # module
