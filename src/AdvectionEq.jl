#===========================================================================================
	Sub-Module		Advection Equation
============================================================================================#
export initialGaussian, rhsAdvection

# The state vector describing how we evolve the advection equation in time
export AdvectionState
struct AdvectionState{T}
    time::T
    u::FinDiff2018.PLCFun{T, T}
end

# State vectors can be scaled and added
function Base. *(a::T, s::AdvectionState{T})::AdvectionState{T} where {T}
    AdvectionState{T}(s.time, a * s.u)
end

function Base. +(s1::AdvectionState{T},
                 s2::AdvectionState{T})::AdvectionState{T} where {T}
    @assert abs(s1.time - s2.time) <= 100*eps(T)
    AdvectionState{T}(s1.time, s1.u + s2.u)
end

############################################################

# A Gaussian, centred at x=1/2, with a width of 1/10
#export guassian
function gaussian(x::T)::T where {T}
    exp(- ((x - 0.5) * 10)^2)
end

# Define initial conditions: A Gaussian at t=0
#export initialGaussian
function initialGaussian(nlines::Int)::AdvectionState{Float64}
    t = 0
    u = samplePLC(Float64, Float64, nlines, gaussian)
    AdvectionState{Float64}(t, u)
end

# The RHS (Right Hand Side) of the advection equation, including boundary conditions
#export rhsAdvection
function rhsAdvection(s::AdvectionState{T})::AdvectionState{T} where {T}
    t = s.time
    u = s.u
    nlines = length(u.points) - 1
    # Calculate spatial derivative
    ux = derivRight(u)
    # Define time derivative as the spatial derivative almost everywhere,
    # except at the right boundary where we set the time derivative to zero
    ut = PLCFun{T, T}([[ux.points[i] for i in 1:nlines]; 0])
    AdvectionState{T}(t, ut)
end

