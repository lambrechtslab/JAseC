__precompile__()

module dataProcessKit
# using SortingAlgorithms
using SortingLab

include("compatible.jl")
export randperm

import Base: float, iterate, length
#Support float("1.23")
float(X::AbstractString) = Base.parse(Float64, X)
float(X::AbstractArray{T}) where {T<:AbstractString} = float.(X)

#Let Regex callable (so filter(r"...", vec) will work).
(R::Regex)(X::AbstractString) = occursin(R, X)

#Let Regex work as a scalar so it can be Broadcast (It's default in Julia 0.6 and also will be a default feature in Julia 1.1.0)
length(::Regex)=1
iterate(R::Regex)=(R, nothing)
iterate(R::Regex, ::Nothing)=nothing
export float, iterate, length

#Support ("String")' = "String"
import Base: adjoint
adjoint(S::Union{AbstractString, Symbol, AbstractChar})=S
export adjoint

includelist=["myEnhanceSyntax.jl"
             "myRowAndGroup.jl"
             "myTable.jl"
             "myCalculation.jl"
             "myBioinformatics.jl"]
for fn in includelist
    include(joinpath(@__DIR__, fn))
end

function iterate(X::Union{eachr, eachcombr, eachgrp, fastgrp}, S=start(X))
    if done(X, S)
        nothing
    else
        next(X, S)
    end
end

# for fn in includelist
#     addhelpfromfile(joinpath(mypth, fn), inmodule=@__MODULE__)
# end

end#module my

