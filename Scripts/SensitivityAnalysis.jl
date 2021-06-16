include("CoinTossXUtilities.jl")
include("ABMVolatilityAuctionProxy.jl")
import HypothesisTests.ADFTest
import Statistics: mean, quantile, std
import StatsBase.kurtosis
import GLM: lm, coef
using ARCHModels, Polynomials, ProgressMeter, CSV#, StatsPlots, DataFrames
#---------------------------------------------------------------------------------------------------

#----- Moments and parameter grid -----#
struct Moments # Moments of log-returns
    μ::Float64 # Mean
    σ::Float64 # Standard deviation
    κ::Float64 # Kurtosis
    hurst::Float64 # Hurst exponent: hurst < 0.5 => mean reverting; hurst == 0.5 => random walk; hurst > 0.5 => momentum
    gph::Float64 # GPH estimator representing long-range dependence
    adf::Float64 # ADF statistic representing random walk property of returns
    garch::Float64 # GARCH paramaters representing short-range dependence
    hill::Float64 # Hill estimator
    function Moments(x)
        logreturns = diff(log.(x))
        μ = mean(logreturns); σ = std(logreturns); κ = kurtosis(logreturns)
        hurst = HurstExponent(logreturns)
        gph = GPH(abs.(logreturns))
        adf = ADFTest(logreturns, :none, 0).stat
        garch = sum(coef(ARCHModels.fit(GARCH{1, 1}, logreturns))[2:3])
        hill = HillEstimator(logreturns[findall(x -> (x >= quantile(logreturns, 0.95)) && (x > 0), logreturns)], 50)
        new(μ, σ, κ, hurst, gph, adf, garch, hill)
    end
end
#---------------------------------------------------------------------------------------------------

#----- Hurst exponent -----#
function HurstExponent(x, d = 50)
    N = length(x)
    if mod(N, 2) != 0 x = push!(x, (x[N - 1] + x[N]) / 2); N += 1 end
    N₁ = N₀ = min(floor(0.99 * N), N-1); dv = Divisors(N₁, d)
    for i in (N₀ + 1):N
        dw = Divisors(i, d)
        if length(dw) > length(dv) N₁ = i; dv = copy(dw) end
    end
    OptN = Int(N₁); d = dv
    x = x[1:OptN]
    RSempirical = map(i -> RS(x, i), d)
    return coeffs(Polynomials.fit(Polynomial, log10.(d), log10.(RSempirical), 1))[2] # Hurst is slope of log-log linear fit
end
function Divisors(n, n₀)
    temp = n₀:floor(n/2)
    return temp[findall(x -> mod(n, x) == 0, temp)]
end
function RS(z, n)
    y = reshape(z, (Int(n), Int(length(z) / n)))
    μ = mean(y, dims = 1)
    σ = std(y, dims = 1)
    temp = cumsum(y .- μ, dims = 1)
    return mean((maximum(temp, dims = 1) - minimum(temp, dims = 1)) / σ)
end
#---------------------------------------------------------------------------------------------------

#----- GPH estimator -----#
function GPH(x, bandwidthExponent = 0.5)
    n = length(x); g = Int(trunc(n^bandwidthExponent))
    j = 1:g; kk = 1:(n - 1)
    w = 2 .* π .* j ./ n # x .-= mean(x)
    σ = sum(x .^ 2) / n
    Σ = map(k -> sum(x[1:(n - k)] .* x[(1 + k):n]) / n, kk)
    periodogram = map(i -> σ + 2 * sum(Σ .* cos.(w[i] .* kk)), j)
    indeces = j[findall(x -> x > 0, periodogram)]
    x_reg = 2 .* log.(2 .* sin.(w[indeces] ./ 2)); y_reg = log.(periodogram[indeces] ./ (2 * π))
    regression = lm(hcat(ones(length(x_reg)), x_reg), y_reg)
    return abs(coef(regression)[2])
end
#---------------------------------------------------------------------------------------------------

#----- Hill estimator -----#
function HillEstimator(x, iterations)
    N = length(x)
    logx = log.(x)
    L = minimum(x); R = maximum(x)
    α = 1 / ((sum(logx) / N) - log(L))
    for i in 1:iterations
        C = (log(L) * (L^(-α)) - log(R) * (R^(-α))) / ((L^(-α)) - (R^(-α)))
        D = (R^α * L^α * (log(L) - log(R))^2) / (L^α - R^α)^2
        α = α * (1 + (α * (sum(logx) / N) - α * C - 1) / (α^2 * D - 1))
    end
    return α
end
#---------------------------------------------------------------------------------------------------

#----- Parameter combinations -----#
Ngrid = [(1, 1, 8), (1, 1, 7), (1, 1, 9), (1, 1, 6), (1, 1, 10)] # 1:4, 2:7, 2:9, 1:3, 1:5
δgrid = range(0.01, 0.1, length = 5)
κgrid = [1, 1.5, 2, 2.5, 3]
νgrid = range(2, 3, length = 5)
σgrid = range(0.1, 1, length = 5)
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

#----- Sensitivity analysis -----#
function SensitivityAnalysis()
    parameterCombinations = CSV.File("/home/ivanjericevich/Parameters.csv") |> DataFrame
    StartJVM()
    gateway = Login(1, 1)
    open("/home/ivanjericevich/Moments.csv", "w") do file
        println(file, "Nt,Nv,Nh,Lambdal,Lambdah,LambdalMin,LambdalMax,LambdahMin,LambdahMax,Delta,Kappa,Nu,M0,Sigma,T,Mu,Sigma,Kappa,Hurst,GPH,ADF,GARCH,Hill") # Header
        @showprogress "Computing:" for i in 3111:nrow(parameterCombinations)
            parameters = Parameters(Nᴸₜ = parameterCombinations[i, :Nt], Nᴸᵥ = parameterCombinations[i, :Nv], Nᴴ = parameterCombinations[i, :Nh], δ = parameterCombinations[i, :Delta], κ = parameterCombinations[i, :Kappa], ν = parameterCombinations[i, :Nu], σ = parameterCombinations[i, :Sigma])
            midPrice = InjectSimulation(gateway, parameters)
            filter!(x -> !ismissing(x), midPrice)
            try
                moments = Moments(midPrice)
                println(file, string(parameters.Nᴸₜ, ",", parameters.Nᴸᵥ, ",", parameters.Nᴴ, ",", parameters.λᴸ, ",", parameters.λᴴ, ",", parameters.λᴸmin, ",", parameters.λᴸmax, ",", parameters.λᴴmin, ",", parameters.λᴴmax, ",", parameters.δ, ",", parameters.κ, ",", parameters.ν, ",", parameters.m₀, ",", parameters.σ, ",", Dates.value(parameters.T), ",", moments.μ, ",", moments.σ, ",", moments.κ, ",", moments.hurst, ",", moments.gph, ",", moments.adf, ",", moments.garch, ",", moments.hill))
            catch e
                println(e)
                println(file, string(parameters.Nᴸₜ, ",", parameters.Nᴸᵥ, ",", parameters.Nᴴ, ",", parameters.λᴸ, ",", parameters.λᴴ, ",", parameters.λᴸmin, ",", parameters.λᴸmax, ",", parameters.λᴴmin, ",", parameters.λᴴmax, ",", parameters.δ, ",", parameters.κ, ",", parameters.ν, ",", parameters.m₀, ",", parameters.σ, ",", Dates.value(parameters.T), ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN))
            end
        end
    end
    Logout(gateway)
    StopCoinTossX()
end
StartCoinTossX(build = false)
SensitivityAnalysis()
#---------------------------------------------------------------------------------------------------

#----- Visualisation -----#
function MomentBoxPlot(; format::String = "pdf")
    data = CSV.File("Data/Moments.csv", types = Dict(:Sigma => String, :Nu => String, :Kappa => String, :Delta => String, :Nh => String)) |> DataFrame |> x -> filter!(y -> !isnan(y.Mu), x)
    for p in [(:Sigma, "σ"), (:Nu, "ν"), (:Kappa, "κ"), (:Delta, "δ"), (:Nh, "N")]
        for m in [(:Mean, "Mean"), (:StandardDeviation, "Standard Deviation"), (:Kurtosis, "Kurtosis"), (:Hurst, "Hurst exponent"), (:GPH, "Geweke and Porter-Hudak estimator"), (:ADF, "
Augmented Dickey–Fuller test statistic"), (:GARCH, "Sum of GARCH(1, 1) parameters"), (:Hill, "Hill estimator")]
            boxplot = @df data violin(p[1], m[1], linewidth = 0, xlabel = p[2], ylabel = m[2])
            @df data boxplot!(boxplot, p[1], m[1], fillalpha = 0.75, linewidth = 2)
            @df data dotplot!(boxplot, p[1], m[1], marker = (:black, stroke(0)))
            savefig(boxplot, string("Figures/", p[1], "-", m[1], "Boxplot.", format))
        end
    end
end
function MomentSurfacePlot(parameter1::Tuple{Symbol, String}, parameter2::Tuple{Symbol, String}, moment::Tuple{Symbol, String}; format::String = "pdf")
    data = CSV.File("Data/Moments.csv") |> DataFrame |> x -> filter!(y -> !isnan(y.Mu), x) |> z -> groupby(z, [parameter1[1], parameter2[1]]) |> a -> combine(a, moment[1] => mean)
    surface = plot(data[:, parameter1[1]], data[:, parameter2[1]], data[:, moment[1]], seriestype = :surface, xlabel = parameter1[2], ylabel = parameter2[2], zlabel = moment[2])
    savefig(surface, string("Figures/", parameter1[1], "-", parameter2[1], "-", moment[1], "Surface.", format))
end
#---------------------------------------------------------------------------------------------------
