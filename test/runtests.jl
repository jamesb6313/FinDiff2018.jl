using Test
using FinDiff2018
#using AdvectionEq
#using TimeEvolce

using Plots

############################################################
##			Piecewise Linear Continuous Functions

# We approximate functions by a piecewise linear continuous approximation, 
# with a regular grid spacing
############################################################

fsinpi = samplePLC(Float64, Float64, 4, sinpi)
#fsinpi2 = samplePLC(Float64, Float64, 4, sinpi)
#@show fsinpi

#@test (fsinpi + fsinpi2) == (fsinpi2 + fsinpi)


(evaluate(fsinpi, 0.3333), sinpi(0.3333))
xs = collect(range(0, stop=1, length=1001))
#@show xs

plot(sinpi.(xs))
plot!([evaluate(fsinpi, x) for x in xs])

rm("C:/Julia/FinDiff2018/test/Plot1_sinPLC.png", force=true)
savefig("Plot1_sinPLC.png")

fsinpi2 = refine(fsinpi)
plot([evaluate(fsinpi, x) for x in xs])
plot!([evaluate(fsinpi2, x)+0.1 for x in xs])

rm("C:/Julia/FinDiff2018/test/Plot1_errorPLC.png", force=true)
savefig("Plot1_errorPLC.png")


############################################################
##			Advection Equation
############################################################
#s0 = initialGaussian(18)
#plot([evaluate(s0.u, x) for x in xs])

#rm("C:/Julia/FinDiff2018/test/Plot2_AdStateInitGuassian.png", force=true)
#savefig("Plot2_AdStateInitGuassian.png")

#r0 = rhsAdvection(s0)
#plot([evaluate(r0.u, x) for x in xs])

#rm("C:/Julia/FinDiff2018/test/Plot2_RHSadvection.png", force=true)
#savefig("Plot2_RHSadvection.png")


############################################################
##			Time evolution
############################################################
#s1 = euler(rhsAdvection, 0.1, s0)

#plot([evaluate(s0.u, x) for x in xs])
#plot!([evaluate(s1.u, x) for x in xs])

#rm("C:/Julia/FinDiff2018/test/Plot3_EulerRHSadvection.png", force=true)
#savefig("Plot3_EulerRHSadvection.png")

#sol1 = solveAdvection(1.0, 100, 0.5);
#sol2 = solveAdvection(1.0, 200, 0.5);

#function plotSolution(sol, every)
#    plot!([evaluate(sol1.states[1].u, x) for x in xs])
#    for state in sol1.states[every:every:end]
#        plot!([evaluate(state.u, x) for x in xs])
#    end
#    plot!()
#end

#plot()
#plotSolution(sol1, 10)

#rm("C:/Julia/FinDiff2018/test/Plot3_solveAdvection.png", force=true)
#savefig("Plot3_solveAdvection.png")

############################################################
##			Discretization Error, Convergence
############################################################


#sol1 = solveAdvection(1.0, 100, 0.5);
#sol2 = solveAdvection(1.0, 200, 0.5);
#sol3 = solveAdvection(1.0, 400, 0.5);

#@show (sol1.states[2].time, sol2.states[3].time, sol3.states[5].time)

#ue1 = refine(sol1.states[2].u) + (-1.0) * sol2.states[3].u
#ue2 = refine(sol2.states[3].u) + (-1.0) * sol3.states[5].u;


#plot([evaluate(ue1, x) for x in xs])
#plot!([evaluate(ue2, x) for x in xs])

#rm("C:/Julia/FinDiff2018/test/Plot4_descreteError_Convergence.png", force=true)
#savefig("Plot4_descreteError_Convergence.png")
