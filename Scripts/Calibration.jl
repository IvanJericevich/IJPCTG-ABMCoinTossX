#=
Calibration:
- Julia version: 1.5.3
- Authors: Ivan Jericevich, Patrick Chang, Tim Gebbie
- Function: Calibrate CoinTossX ABM simulated mid-prices to JSE mid-price data
- Structure:
    1. Moving block bootstrap to estimate covariance matrix of empirical moments on JSE mid-price time-series
    2. Objective function to be minimized
    3. Calibrate with NMTA optimization
    4. Visualisation
- Examples:
    midprice = CSV.File("Data/JSECleanedTAQNPN.csv", select = [:MidPrice], limit = 40000) |> Tables.matrix |> vec |> x -> filter(y -> !isnan(y), x)
    W = MovingBlockBootstrap(midprice, 500)
    Calibrate(initialsolution)
=#
using JLD, CSV, Plots
import Statistics: cov
addprocs(2)
@everywhere import Random.rand
@everywhere include("Scripts/Moments.jl")
include("ABMVolatilityAuctionProxy.jl"); include("NMTA.jl")
#---------------------------------------------------------------------------------------------------

#----- Moving block bootstrap to estimate covariance matrix of empirical moments on JSE mid-price time-series -----#
function MovingBlockBootstrap(midprice::Vector{Float64}, iterations::Int64 = 1000, windowsize::Int64 = 2000)
    logreturns = diff(log.(midprice))
    bootstrapmoments = @distributed (vcat) for i in 1:iterations
        indeces = rand(1:(length(logreturns) - windowsize + 1), Int(ceil(length(logreturns)/windowsize)))
        bootstrapreturns = Vector{Float64}()
        for index in indeces
            append!(bootstrapreturns, logreturns[index:(index  + windowsize - 1)])
        end
        moments = Moments(bootstrapreturns[1:length(logreturns)], logreturns)
        [moments.μ moments.σ moments.κ moments.ks moments.hurst moments.gph moments.adf moments.garch moments.hill]
    end
    W = inv(cov(bootstrapmoments))
    save("Data/W.jld", "W", W)
end
#---------------------------------------------------------------------------------------------------

#----- Objective function to be minimized -----#
function WeightedSumofSquaredErrors(parameters::Parameters, replications::Int64, W::Array{Float64, 2}, empiricalmoments::Moments, empiricallogreturns::Vector{Float64}, gateway::TradingGateway)
    errormatrix = fill(0.0, (replications, 9))
    for i in 1:replications
        microprice = InjectSimulation(gateway, parameters, seed = i)
        if !isempty(microprice)
            filter!(x -> !isnan(x), microprice)
            logreturns = diff(log.(microprice))
            try
                simulatedmoments = Moments(logreturns, empiricallogreturns)
                errormatrix[i, :] = [simulatedmoments.μ-empiricalmoments.μ simulatedmoments.σ-empiricalmoments.σ simulatedmoments.κ-empiricalmoments.κ simulatedmoments.ks-empiricalmoments.ks simulatedmoments.hurst-empiricalmoments.hurst simulatedmoments.gph-empiricalmoments.gph simulatedmoments.adf-empiricalmoments.adf simulatedmoments.garch-empiricalmoments.garch simulatedmoments.hill-empiricalmoments.hill]
            catch e
                println(e)
                errormatrix[i, :] = errormatrix[i - 1, :]
            end
        else
            return Inf
        end
    end
    GC.gc() # Garbage collection
    errors = mean(errormatrix, dims = 1)
    return (errors * W * transpose(errors))[1]
end
#---------------------------------------------------------------------------------------------------

#----- Calibrate with NMTA optimization -----#
function Calibrate(initialsolution::Vector{Float64}; f_reltol::Vector{Float64} = [0.3, 0.2, 0.1, 0], ta_rounds::Vector{Int64} = [12, 10, 8, 6], neldermeadstate = nothing)
    empiricallogreturns = CSV.File("JSE/L1LOB.csv", missingstring = "missing", ignoreemptylines = true, select = [:MicroPrice], skipto = 20000, limit = 20000) |> Tables.matrix |> vec |> y -> filter(z -> !ismissing(z), y) |> x -> diff(log.(x))
    empiricalmoments = Moments(empiricallogreturns, empiricallogreturns)
    StartCoinTossX(build = false); StartJVM(); sleep(20); gateway = Login(1, 1)
    W = load("W.jld")["W"]
    counter = Counter(0)
    objective = NonDifferentiable(x -> WeightedSumofSquaredErrors(Parameters(Nᴴ = Int(ceil(x[1])), δ = abs(x[2]), κ = abs(x[3]), ν = abs(x[4]), σ = abs(x[5]), T = Millisecond(3600 * 2 * 1000)), 4, W, empiricalmoments, empiricallogreturns, gateway), initialsolution)
    optimizationoptions = Options(show_trace = true, store_trace = true, trace_simplex = true, extended_trace = true, iterations = !isempty(ta_rounds) ? sum(ta_rounds) : 30, ξ = 0.15, ta_rounds = ta_rounds, f_reltol = f_reltol)
    result = !isnothing(neldermeadstate) ? Optimize(objective, initialsolution, optimizationoptions, neldermeadstate) : Optimize(objective, initialsolution, optimizationoptions)
    save("OptimizationResult.jld", "result", result)
    Logout(gateway); StopCoinTossX()
end
#---------------------------------------------------------------------------------------------------

#----- Validate optimization results -----#
stacktrace = load("Data/Calibration/OptimizationResults.jld")["result"]
f = zeros(Float64, length(stacktrace)); g_norm = zeros(Float64, length(stacktrace)); f_simplex = fill(0.0, length(stacktrace), 6)#; centr = fill(0.0, length(stacktrace), 5); metadata = Vector{Dict{Any, Any}}()
for i in 1:length(stacktrace)
    f[i] = stacktrace[i].value
    g_norm[i] = stacktrace[i].g_norm
    push!(metadata, stacktrace[i].metadata)
    f_simplex[i, :] = transpose(stacktrace[i].metadata["simplex_values"])
end
# Objectives
objectives = plot(1:length(stacktrace), f, seriestype = :line, linecolor = :blue, label = "Weighted SSE objective", xlabel = "Iteration", ylabel = "Weighted SSE objective", legendfontsize = 5, fg_legend = :transparent, tickfontsize = 5, xaxis = false, xticks = false, legend = :bottomleft, guidefontsize = 7, yscale = :log10, minorticks = true)
plot!(twinx(), 1:length(stacktrace), g_norm, seriestype = :line, linecolor = :purple, label = "Convergence criterion", ylabel = "Convergence criterion", legend = :topright, legendfontsize = 5, fg_legend = :transparent, tickfontsize = 5, yscale = :log10, minorticks = true)
savefig(objectives, "Figures/NMTAFitness.pdf")
# Simplex values
convergence = plot(1:length(stacktrace), f_simplex, seriestype = :line, linecolor = [:blue :purple :green :orange :red :black], xlabel = "Iteration", ylabel = "Weighted SSE objective", legend = false, tickfontsize = 5, guidefontsize = 7, yscale = :log10, minorticks = true)
savefig(convergence, "Figures/NMTASimplexFitnesses.pdf")
#---------------------------------------------------------------------------------------------------
