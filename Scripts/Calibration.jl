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
    Calibrate(W, midprice)
=#
using DataFrames, JLD, Plots#, Distributed
import CSV: File
import Random.rand
# import Statistics: cov
# addprocs(2)
# @everywhere import Random.rand
# @everywhere include("Scripts/Moments.jl")
include("CoinTossXUtilities.jl"); include("ABMVolatilityAuctionProxy.jl"); include("NMTA.jl")
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
            simulatedmoments = Moments(logreturns, empiricallogreturns)
            errormatrix[i, :] = [simulatedmoments.μ-empiricalmoments.μ simulatedmoments.σ-empiricalmoments.σ simulatedmoments.κ-empiricalmoments.κ simulatedmoments.ks-empiricalmoments.ks simulatedmoments.hurst-empiricalmoments.hurst simulatedmoments.gph-empiricalmoments.gph simulatedmoments.adf-empiricalmoments.adf simulatedmoments.garch-empiricalmoments.garch simulatedmoments.hill-empiricalmoments.hill]
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
function Calibrate(initialsolution::Vector{Float64}, f_reltol::Vector{Float64} = [0.3, 0.2, 0.1, 0], ta_rounds::Vector{Int64} = [14, 12, 10, 8], neldermeadstate = nothing)
    empiricallogreturns = File("Data/JSEL1LOB.csv", missingstring = "missing", ignoreemptylines = true, select = [:MicroPrice], skipto = 20000, limit = 20000) |> Tables.matrix |> vec |> y -> filter(z -> !ismissing(z), y) |> x -> diff(log.(x))
    empiricalmoments = Moments(empiricallogreturns, empiricallogreturns)
    StartCoinTossX(build = false); StartJVM(); sleep(20); gateway = Login(1, 1)
    W = load("Data/W.jld")["W"]

    try
        objective = NonDifferentiable(x -> WeightedSumofSquaredErrors(Parameters(Nᴴ = Int(ceil(x[1])), δ = x[2], κ = x[3], ν = x[4], σ = x[5], T = Millisecond(3600 * 2 * 1000)), 10, W, empiricalmoments, empiricallogreturns, gateway), initialsolution)
        optimizationoptions = Options(show_trace = true, store_trace = true, trace_simplex = true, extended_trace = true, iterations = sum(ta_rounds), ξ = 0.15, ta_rounds = ta_rounds, f_reltol = f_reltol)
        result = !isnothing(neldermeadstate) ? Optimize(objective, initialsolution, optimizationoptions, neldermeadstate) : Optimize(objective, initialsolution, optimizationoptions)
        save("Data/OptimizationResults.jld", "result", result)
    catch e
        println(e)
        Logout(gateway); StopCoinTossX()
    end
    Logout(gateway); StopCoinTossX()
end
initialsolution = load("BestParameters.jld")["parameters"][end - 1]
neldermeadstate = end_state(load("OptimizationResults.jld")["result"])
Calibrate([initialsolution.Nᴴ, initialsolution.δ, initialsolution.κ, initialsolution.ν, initialsolution.σ], 10, neldermeadstate)
#---------------------------------------------------------------------------------------------------

#----- Validate optimization results -----#
result = load("Data/OptimizationResults.jld")["result"]
typeof(result)
stacktrace = trace(result)
f = zeros(Float64, length(trace)); g_norm = zeros(Float64, length(trace)); f_simplex = fill(0.0, length(trace), 3); centroid = fill(0.0, length(trace), 2)
for i in 1:length(trace)
    f[i] = trace[i].value
    g_norm[i] = trace[i].g_norm
    metadata = trace[i].metadata
    f_simplex[i, :] = transpose(metadata["simplex_values"])
    centroid[i, :] = transpose(metadata["centroid"])
end
# Objectives
objectives = plot(1:length(trace), f, seriestype = :line, linecolor = :blue, label = "Weighted SSE objective", xlabel = "Iteration", ylabel = "Value", legendfontsize = 5, fg_legend = :transparent, tickfontsize = 5, xaxis = false, xticks = false)
plot!(twinx(), 1:length(trace), g_norm, seriestype = :line, linecolor = :purple, label = "Convergence criterion", legend = :topright)
# Simplex values
plot(1:length(trace), f_simplex, seriestype = :line, linecolor = [:blue :purple :green], xlabel = "Iteration", ylabel = "Simplex objective values")
# Centroids
plot(1:length(trace), centroid, seriestype = :line, linecolor = [:blue :purple], xlabel = "Iteration", ylabel = "Simplex objective values")
#---------------------------------------------------------------------------------------------------
