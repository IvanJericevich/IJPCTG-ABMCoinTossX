#=
SensitivityAnalysis:
- Julia version: 1.5.3
- Authors: Ivan Jericevich, Patrick Chang, Tim Gebbie
- Function: Simulated grid of moments for all combinations of realistic parameters
- Structure:
    1. Sensitivity analysis
    2. Parameter combinations
    3. Visualisation
- Examples:
    StartCoinTossX(build = false)
    SensitivityAnalysis()
=#
using ProgressMeter, CSV, JLD, Plots, DataFrames, StatsPlots, Statistics, ColorSchemes, Dates
using LinearAlgebra: diag, inv
include("ABMVolatilityAuctionProxy.jl"); include("Moments.jl")
#---------------------------------------------------------------------------------------------------

#----- Sensitivity analysis -----#
function SensitivityAnalysis()
    parameterCombinations = CSV.File("Data/Sensitivity/Parameters.csv") |> DataFrame # Ignore header
    empiricallogreturns = load("Data/EmpiricalLogReturns.jld")["empiricallogreturns"]
    empiricalmoments = load("Data/EmpiricalMoments.jld")["empiricalmoments"]
    W = load("Data/W.jld")["W"]
    StartJVM()
    gateway = Login(1, 1)
    open(string("Results.csv"), "w") do file
        println(file, "Side,Nt,Nv,Nh,Lambdal,Lambdah,LambdalMin,LambdalMax,LambdahMin,LambdahMax,Delta,Kappa,Nu,M0,Sigma,T,Mean,StandardDeviation,Kurtosis,KolmogorovSmirnov,Hurst,GPH,ADF,GARCH,Hill,Objective,RunTime") # Header
        @showprogress "Computing:" for i in 1:nrow(parameterCombinations)
            parameters = Parameters(Nᴸₜ = parameterCombinations[i, :Nt], Nᴸᵥ = parameterCombinations[i, :Nv], Nᴴ = parameterCombinations[i, :Nh], δ = parameterCombinations[i, :Delta], κ = parameterCombinations[i, :Kappa], ν = parameterCombinations[i, :Nu], σ = parameterCombinations[i, :Sigma], T = Millisecond(3600 * 1000))
            t₀ = time()
            try
                midPrice, microPrice = InjectSimulation(gateway, parameters)
                filter!(x -> !isnan(x), midPrice); filter!(x -> !isnan(x), microPrice)
                midpricemoments = Moments(diff(log.(midPrice)), first(empiricallogreturns)); micropricemoments = Moments(diff(log.(microPrice)), last(empiricallogreturns))
                midpriceerror = [midpricemoments.μ-first(empiricalmoments).μ, midpricemoments.σ-first(empiricalmoments).σ, midpricemoments.κ-first(empiricalmoments).κ, midpricemoments.ks-first(empiricalmoments).ks, midpricemoments.hurst-first(empiricalmoments).hurst, midpricemoments.gph-first(empiricalmoments).gph, midpricemoments.adf-first(empiricalmoments).adf, midpricemoments.garch-first(empiricalmoments).garch, midpricemoments.hill-first(empiricalmoments).hill]
                micropriceerror = [micropricemoments.μ-last(empiricalmoments).μ, micropricemoments.σ-last(empiricalmoments).σ, micropricemoments.κ-last(empiricalmoments).κ, micropricemoments.ks-last(empiricalmoments).ks, micropricemoments.hurst-last(empiricalmoments).hurst, micropricemoments.gph-last(empiricalmoments).gph, micropricemoments.adf-last(empiricalmoments).adf, micropricemoments.garch-last(empiricalmoments).garch, micropricemoments.hill-last(empiricalmoments).hill]
                midpriceobjective = transpose(midpriceerror) * W * midpriceerror; micropriceobjective = transpose(micropriceerror) * W * micropriceerror
                println(file, string("MidPrice,", parameters.Nᴸₜ, ",", parameters.Nᴸᵥ, ",", parameters.Nᴴ, ",", parameters.λᴸ, ",", parameters.λᴴ, ",", parameters.λᴸmin, ",", parameters.λᴸmax, ",", parameters.λᴴmin, ",", parameters.λᴴmax, ",", parameters.δ, ",", parameters.κ, ",", parameters.ν, ",", parameters.m₀, ",", parameters.σ, ",", Dates.value(parameters.T), ",", midpricemoments.μ, ",", midpricemoments.σ, ",", midpricemoments.κ, ",", midpricemoments.ks, ",", midpricemoments.hurst, ",", midpricemoments.gph, ",", midpricemoments.adf, ",", midpricemoments.garch, ",", midpricemoments.hill, ",", midpriceobjective, ",", time() - t₀))
                println(file, string("MicroPrice,", parameters.Nᴸₜ, ",", parameters.Nᴸᵥ, ",", parameters.Nᴴ, ",", parameters.λᴸ, ",", parameters.λᴴ, ",", parameters.λᴸmin, ",", parameters.λᴸmax, ",", parameters.λᴴmin, ",", parameters.λᴴmax, ",", parameters.δ, ",", parameters.κ, ",", parameters.ν, ",", parameters.m₀, ",", parameters.σ, ",", Dates.value(parameters.T), ",", micropricemoments.μ, ",", micropricemoments.σ, ",", micropricemoments.κ, ",", micropricemoments.ks, ",", micropricemoments.hurst, ",", micropricemoments.gph, ",", micropricemoments.adf, ",", micropricemoments.garch, ",", micropricemoments.hill, ",", micropriceobjective, ",", time() - t₀))
            catch e
                println(e)
                println(file, string("MidPrice,", parameters.Nᴸₜ, ",", parameters.Nᴸᵥ, ",", parameters.Nᴴ, ",", parameters.λᴸ, ",", parameters.λᴴ, ",", parameters.λᴸmin, ",", parameters.λᴸmax, ",", parameters.λᴴmin, ",", parameters.λᴴmax, ",", parameters.δ, ",", parameters.κ, ",", parameters.ν, ",", parameters.m₀, ",", parameters.σ, ",", Dates.value(parameters.T), ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN))
                println(file, string("MicroPrice,", parameters.Nᴸₜ, ",", parameters.Nᴸᵥ, ",", parameters.Nᴴ, ",", parameters.λᴸ, ",", parameters.λᴴ, ",", parameters.λᴸmin, ",", parameters.λᴸmax, ",", parameters.λᴴmin, ",", parameters.λᴴmax, ",", parameters.δ, ",", parameters.κ, ",", parameters.ν, ",", parameters.m₀, ",", parameters.σ, ",", Dates.value(parameters.T), ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN))
            end
            GC.gc()
        end
    end
    Logout(gateway)
    StopCoinTossX()
end
#---------------------------------------------------------------------------------------------------

#----- Parameter combinations -----#
Ngrid = [(1, 1, 14), (1, 1, 13), (1, 1, 12), (1, 1, 11), (1, 1, 10)]
δgrid = range(0.01, 0.17, length = 5)
κgrid = [1.5, 2, 2.5, 3, 3.5]
νgrid = range(2, 4.8, length = 5)
σgrid = range(0.01, 1, length = 5)
open("/home/ivanjericevich/Parameters.csv", "w") do file
    println(file, "Nt,Nv,Nh,Lambdal,Lambdah,LambdalMin,LambdalMax,LambdahMin,LambdahMax,Delta,Kappa,Nu,M0,Sigma,T") # Header
    for N in Ngrid
        for δ in δgrid
            for κ in κgrid
                for ν in νgrid
                    for σ in σgrid
                        parameters = Parameters(Nᴸₜ = N[1], Nᴸᵥ = N[2], Nᴴ = N[3], δ = δ, κ = κ, ν = ν, σ = σ)
                        println(file, string(parameters.Nᴸₜ, ",", parameters.Nᴸᵥ, ",", parameters.Nᴴ, ",", parameters.λᴸ, ",", parameters.λᴴ, ",", parameters.λᴸmin, ",", parameters.λᴸmax, ",", parameters.λᴴmin, ",", parameters.λᴴmax, ",", parameters.δ, ",", parameters.κ, ",", parameters.ν, ",", parameters.m₀, ",", parameters.σ, ",", Dates.value(parameters.T)))
                    end
                end
            end
        end
    end
end
#---------------------------------------------------------------------------------------------------

#----- Visualisation -----#
function MomentBoxPlot(timeseries::Symbol; format::String = "pdf")
    data = CSV.File("Data/Sensitivity/Moments.csv", types = Dict(:Series => Symbol)) |> DataFrame |> x -> filter!(y -> !ismissing(y.Objective) && y.Series == timeseries, x)
    parameters = [(:Sigma, "σ"), (:Nu, "ν"), (:Kappa, "κ"), (:Delta, "δ"), (:Nh, "N")]
    moments = [(:Kurtosis, "Kurtosis"), (:Hurst, "Hurst exponent"), (:GPH, "GPH estimator"), (:ADF, "ADF test statistic"), (:Hill, "Hill estimator")]
    matrixPlot = plot(layout = (length(parameters), length(moments)), xrotation = 30, yrotation = 30, guidefontsize = 5, tickfontsize = 3, tick_direction = :out, ylabel = ["Kurtosis" "" "" "" "" "Hurst exponent" "" "" "" "" "GPH estimator" "" "" "" "" "ADF statistic" "" "" "" "" "Hill estimator" "" "" "" ""], xlabel = ["" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "σ" "ν" "κ" "δ" "N"])
    counter = 1
    for m in 1:length(moments)
        for p in 1:length(parameters)
            violin!(matrixPlot, data[:, parameters[p][1]], data[:, moments[m][1]], linealpha = 0, fillalpha = 0.8, fillcolor = :blue, subplot = counter, legend = false)
            boxplot!(matrixPlot, data[:, parameters[p][1]], data[:, moments[m][1]], fillalpha = 0, marker = (1, :black, stroke(:black)), linewidth = 0, linecolor = :black, subplot = counter, legend = false)
            counter += 1
        end
    end
    savefig(matrixPlot, string("Figures/", timeseries, "MatrixBoxplot.", format))
end
function MomentSurface(timeseries::Symbol, moment::Tuple{Symbol, String}; format::String = "pdf")
    parameters = [(:Sigma, "σ"), (:Nu, "ν"), (:Kappa, "κ"), (:Delta, "δ"), (:Nh, "N")]
    data = CSV.File("Data/Sensitivity/Moments.csv", types = Dict(:Series => Symbol)) |> DataFrame |> x -> filter!(y -> !ismissing(y.Objective) && y.Series == timeseries, x)
    surface1 = groupby(data, first.(parameters[[4, 3]])) |> a -> combine(a, moment[1] => mean) |> df -> plot(unique(df[:, parameters[4][1]]), unique(df[:, parameters[3][1]]), transpose(reshape(df[:, string(moment[1], "_mean")], (5, 5))), seriestype = :surface, xlabel = parameters[4][2], ylabel = parameters[3][2], zlabel = moment[2], colorbar = false)
    surface2 = groupby(data, first.(parameters[[5, 2]])) |> a -> combine(a, moment[1] => mean) |> df -> plot(unique(df[:, parameters[5][1]]), unique(df[:, parameters[2][1]]), transpose(reshape(df[:, string(moment[1], "_mean")], (5, 5))), seriestype = :surface, xlabel = parameters[5][2], ylabel = parameters[2][2], zlabel = moment[2], colorbar = false)
    surface = plot(surface1, surface2, layout = (1, 2))
    savefig(surface, string("Figures/", timeseries, moment[1], "Surface.", format))
end
function ObjectiveSurface(timeseries::Symbol; format::String = "pdf")
    parameters = [(:Sigma, "σ"), (:Nu, "ν"), (:Kappa, "κ"), (:Delta, "δ"), (:Nh, "N")]
    data = CSV.File("Data/Sensitivity/Moments.csv", select = vcat([:Series, :Objective], first.(parameters)), types = Dict(:Series => Symbol)) |> DataFrame |> x -> filter!(y -> !ismissing(y.Objective) && y.Series == timeseries, x)
    surface1 = groupby(data, first.(parameters[[4, 3]])) |> a -> combine(a, :Objective => mean) |> df -> plot(unique(df[:, first(parameters[4])]), unique(df[:, first(parameters[3])]), transpose(reshape(df[:, "Objective_mean"], (5, 5))), seriestype = :surface, xlabel = last(parameters[4]), ylabel = last(parameters[3]), zlabel = "Objective", colorbar = false, fill = :seismic)
    surface2 = groupby(data, first.(parameters[[5, 2]])) |> a -> combine(a, :Objective => mean) |> df -> plot(unique(df[:, first(parameters[5])]), unique(df[:, first(parameters[2])]), transpose(reshape(df[:, "Objective_mean"], (5, 5))), seriestype = :surface, xlabel = last(parameters[5]), ylabel = last(parameters[2]), zlabel = "Objective", colorbar = false, fill = :seismic)
    surface = plot(surface1, surface2, layout = (1, 2))
    savefig(surface, string("Figures/", timeseries, "ObjectiveSurface.", format))
end
function ConfidenceIntervals(; format::String = "pdf")
    variables = [(:Nh, "N"), (:Delta, "δ"), (:Kappa, "κ"), (:Nu, "ν"), (:Sigma, "σ"), (:Mean, "Mean"), (:StandardDeviation, "StDev."), (:Kurtosis, "Kurtosis"), (:KolmogorovSmirnov, "KS"), (:Hurst, "Hurst"), (:GPH, "GPH"), (:ADF, "ADF"), (:GARCH, "GARCH"), (:Hill, "Hill")]
    W = load("Data/Calibration/W.jld")["W"]
    # Parameters
    data = CSV.File("Data/Sensitivity/Moments.csv", select = vcat([:Series], first.(variables)), types = Dict(:Series => Symbol)) |> DataFrame |> x -> filter!(y -> !ismissing(y.Mean) && y.Series == :MicroPrice, x) |> df -> Matrix(df[:, 2:end])
    β = cov(data)[1:5, 6:end] ./ var(data[:, 6:end], dims = 1)
    Σ = β * inv(W) * transpose(β)
    σ = sqrt.(diag(Σ))
    estimates = minimizer(load("Data/Calibration/OptimizationResult.jld")["result"])
    upper = estimates .+ (1.96 .* σ)
    lower = estimates .- (1.96 .* σ)
    open("Data/Parameters.txt", "w") do file
        println(file, "Variable,Lower,Estimate,Upper")
        for i in 1:5
            if i == 1
                println(file, string(first(variables[i]), ",", round(lower[i] / 2, digits = 3), ",", round(estimates[i] / 2, digits = 3), ",", round(upper[i] / 2, digits = 3)))
            else
                println(file, string(first(variables[i]), ",", round(lower[i], digits = 3), ",", round(estimates[i], digits = 3), ",", round(upper[i], digits = 3)))
            end
        end
    end
    Ρ = cor(data)
    correlations = heatmap(last.(variables), last.(variables), Ρ, c = cgrad(:seismic, [0, 0.25, 0.45, 1]), xrotation = 30, yrotation = 30)
    annotate!(correlations, [(j - 0.5, i - 0.5, text(round(Ρ[i,j], digits = 3), 5, :black, :center)) for i in 1:size(Ρ, 1) for j in 1:size(Ρ, 1)], linecolor = :white)
    savefig(correlations, string("Figures/ParameterMomentCorrelation.", format))
    weights = heatmap(last.(variables[6:end]), last.(variables[6:end]), log.(W), c = cgrad(:seismic, [0, 0.25, 0.45, 1]), xrotation = 30, yrotation = 30)
    annotate!(weights, [(j - 0.5, i - 0.5, text(round(log(W[i,j]), digits = 3), 5, :black, :center)) for i in 1:size(W, 1) for j in 1:size(W, 1)], linecolor = :white)
    savefig(weights, string("Figures/InverseCovarianceEmpiricalMoments.", format))
    # Moments
    cointossx = CSV.File(string("Data/Calibration/L1LOB.csv"), select = [:DateTime, :MicroPrice], missingstring = "missing") |> DataFrame
    jse = CSV.File(string("Data/JSE/L1LOB.csv"), select = [:DateTime, :MicroPrice], missingstring = "missing") |> DataFrame |> x -> filter!(y -> y.DateTime < DateTime("2019-01-04T09:00:41"), x)
    cointossx.Date = Date.(cointossx.DateTime); jse.Date = Date.(jse.DateTime)
    cointossxUniqueDays = unique(cointossx.Date); jseUniqueDays = unique(jse.Date)
    cointossxlogreturns = map(day -> diff(log.(skipmissing(cointossx[searchsorted(cointossx.Date, day), :MicroPrice]))), cointossxUniqueDays) |> x -> reduce(vcat, x)
    jselogreturns = map(day -> diff(log.(skipmissing(jse[searchsorted(jse.Date, day), :MicroPrice]))), jseUniqueDays) |> x -> reduce(vcat, x)
    cointossxmoments = Moments(cointossxlogreturns, jselogreturns); jsemoments = Moments(jselogreturns, jselogreturns)
    empirical = [jsemoments.μ, jsemoments.σ, jsemoments.κ, jsemoments.ks, jsemoments.hurst, jsemoments.gph, jsemoments.adf, jsemoments.garch, jsemoments.hill]
    estimates = [cointossxmoments.μ, cointossxmoments.σ, cointossxmoments.κ, cointossxmoments.ks, cointossxmoments.hurst, cointossxmoments.gph, cointossxmoments.adf, cointossxmoments.garch, cointossxmoments.hill]
    σ = sqrt.(diag(inv(W)))
    upper = estimates .+ (1.96 .* σ)
    lower = estimates .- (1.96 .* σ)
    open("Data/Moments.txt", "w") do file
        println(file, "Moment,Lower,Estimate,Upper,Empirical")
        for i in 1:length(estimates)
            println(file, string(first(variables[i+5]), ",", round(lower[i], digits = 5), ",", round(estimates[i], digits = 5), ",", round(upper[i], digits = 5), ",", round(empirical[i], digits = 5)))
        end
    end
end
#---------------------------------------------------------------------------------------------------
