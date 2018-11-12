module FinDiff2018

#push!(LOAD_PATH, @__DIR__) # expose all other modules defined in this directory.

#using Plots

export PLCFun, derivRight, samplePLC, evaluate, refine

############################################################
##			Piecewise Linear Continuous Functions

# We approximate functions by a piecewise linear continuous approximation, 
# with a regular grid spacing
############################################################


# A piecewise linear continuous function mapping a type T to a type U
#export PLCFun
struct PLCFun{T,U}
    # Domain: [0; 1]
    points::Vector{U}
end

# Functions can be scaled and added
# function Base.length(f::PLCFun{T, U})::Int where {T, U}
#     length(f.points)
# end
function Base. *(a::U, f::PLCFun{T, U})::PLCFun{T, U} where {T, U}
    PLCFun{T, U}(a .* f.points)
end
function Base. +(f::PLCFun{T, U}, g::PLCFun{T, U})::PLCFun{T, U} where {T, U}
    @assert length(f.points) == length(g.points)
    PLCFun{T, U}(f.points .+ g.points)
end


############################################################
# Calculate the x coordinates of the endpoints of the lines
#export xcoord
function xcoord(::Type{T}, nlines::Int, i::Int)::T where {T}
    @assert 1 <= i <= nlines + 1
    dx = 1 / nlines
    x = (i-1) * dx
    x
end


# The inverse of "xcoord": Determine the line segment on which a particular x coordinate lies
#export lineidx
function lineidx(f::PLCFun{T, U}, x::T)::Int where {T, U}
    @assert 0 <= x <= 1
    nlines = length(f.points) - 1
    dx = 1 / nlines
    i = floor(Int, x / dx) + 1
    i = max(1, i)
    i = min(nlines, i)
    i
end


# Convert a general Julia function into a PLCFun. We need to specify the types T and U,
# as well as the number of line segments to use.
#export samplePLC
function samplePLC(::Type{T}, ::Type{U}, nlines::Int, f::Function)::PLCFun{T, U} where {T, U}
    ys = U[f(xcoord(T, nlines, i)) for i in 1:nlines+1]
    PLCFun{T, U}(ys)
end


# Linear interpolation between two points
#export linterp
function linterp(x1::T, y1::U, x2::T, y2::U, x::T)::U where {T, U}
    y1 * (x - x2) / (x1 - x2) + y2 * (x - x1) / (x2 - x1)
end


# Evaluate a PLCFun at point x
#export evaluate
function evaluate(f::PLCFun{T,U}, x::T)::U where {T, U}
    @assert 0 <= x <= 1
    nlines = length(f.points) - 1
    # Find out on which line segment the point lies
    i = lineidx(f, x)
    # Interpolate between the two endpoints of the line segment
    x1 = xcoord(T, nlines, i)
    x2 = xcoord(T, nlines, i+1)
    y1 = f.points[i]
    y2 = f.points[i+1]
    y = linterp(x1, y1, x2, y2, x)
    y
end

# Calculate the derivative of a PLCFun, using a right-biased derivative
#export derivRight
function derivRight(f::PLCFun{T, U})::PLCFun{T, U} where {T, U}
    nlines = length(f.points) - 1
    dx = 1 / nlines
    ys = [[(f.points[i+1] - f.points[i]) / dx for i in 1:nlines];
          (f.points[end] - f.points[end-1]) / dx]
    PLCFun{T, U}(ys)
end

#export refine
function refine(f::PLCFun{T, U})::PLCFun{T, U} where {T, U}
    nlines = length(f.points) - 1
    dx = 1 / nlines
    ys = [evaluate(f, xcoord(T, 2*nlines, i)) for i in 1:2*nlines+1]
    PLCFun{T, U}(ys)
end

#=

Create Plots
#export testplot
#function testplot(arr::Array{Float64,2}, filename::String)
#	heatmap(arr, clim=(-1.0, +1.0), color=:viridis)
#	savefig(filename)
#end

=#


include("AdvectionEq.jl")
include("TimeEvolve.jl")

end # module FinDiff2018