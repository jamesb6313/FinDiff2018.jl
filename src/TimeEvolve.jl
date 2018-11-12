#============================================================================================
			Time Evolve
============================================================================================#

export euler, Solution, solveAdvection

#export euler
function euler(rhs::Function, dt::T,
               s0::AdvectionState{T})::AdvectionState{T} where {T}
    r0 = rhs(s0)
    s1 = s0 + dt * r0
    AdvectionState{T}(s0.time + dt, s1.u)
end


struct Solution{T}
    dx::T
    dt::T
    states::Vector{AdvectionState{T}}
end

function solveAdvection(tmax::T, nlines::Int, lambda::T)::Solution{T} where {T}
    dx = 1 / nlines
    dt = lambda * dx
    nsteps = round(Int, tmax / dt)
    sol = Solution{T}(dx, dt, AdvectionState{T}[])
    s = initialGaussian(nlines)
    push!(sol.states, s)
    for step in 1:nsteps
        s = euler(rhsAdvection, dt, s)
        push!(sol.states, s)
    end
    sol
end