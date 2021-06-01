using HypothesisTests, Statistics, StatsBase, Polynomials, GARCH
include("CoinTossXUtilities.jl")
include("ABM.jl")
#---------------------------------------------------------------------------------------------------

#----- Moments and parameter grid -----#
struct Moments # Moments of log-returns
    μ::Float64 # Mean
    σ::Float64 # Standard deviation
    κ::Float64 # Kurtosis
    α::Float64 # Power-law parameter calibrated to ACF
    hurst::Float64 # Hurst exponent: hurst < 0.5 => mean reverting; hurst == 0.5 => random walk; hurst > 0.5 => momentum
    gph::Float64 # GPH estimator representing long-range dependence
    adf::Float64 # ADF statistic representing random walk property of returns
    garch::Float64 # GARCH paramaters representing short-range dependence
    hill::Float64 # Hill estimator
    function Moments(x::Vector{Float64})
        logreturns = diff(log.(x))
        μ = mean(logreturns); σ = std(logreturns); κ = kurtosis(logreturns)
        α =
        hurst = HurstExponent(midPrice)
        gph =
        adf = ADFTest(logreturns, deterministic, lag)
        garch = sum(garchFit(logreturns).params[2:3])
        hill = HillEstimator(logreturns, iterations)
        new(μ, σ, κ, α, hurst, gph, adf, garch, hill)
    end
end
δgrid
κgrid
νgrid
σgrid
#---------------------------------------------------------------------------------------------------

#----- δ cube -----#
δmoments = Dict{Float64, Array{Float64, 3}}() # Each key is a single value of δ
for δ in δgrid
    combinations = zeros(Float64, length(κgrid), length(νgrid), length(σgrid))
    for x in 1:length(κgrid)
        for y in 1:length(νgrid)
            for z in 1:length(σgrid)
                parameters = Parameters()
                time, midPrice = InjectSimulation()
                combinations[x, y, z] = Moments(midPrice)
            end
        end
    end
    δmoments[δ] = combinations
end
#---------------------------------------------------------------------------------------------------

#----- Hurst exponent -----#
function HurstExponent(x, d = 50)
    N = length(x)
    if mod(N, 2) != 0 x = push!(x, (x[N - 1] + x[N]) / 2); N += 1 end
    N₁ = N₀ = min(floor(0.99 * N), N-1); dv = Divisors(N₁, d)
    for i in (N₀ + 1):N
        dw = Divisors(i, d)
        if length(dw) > length(dv) N₁ = i; dv = dw end
    end
    OptN = Int(N₁); d = dv
    x = x[1:OptN]
    RSempirical = map(i -> RS(x, i), d)
    return coeffs(fit(Polynomial, log10.(d), log10.(RSempirical), 1))[1]
end
function Divisors(n, n₀)
    temp = n₀:floor(n/2)
    return temp[findall(x -> mod(n, x) == 0, temp)]
end
function RS(z, n)
    y = reshape(x, (Int(n), Int(length(z) / n)))
    μ = mean(y, dims = 1)
    σ = std(y, dims = 1)
    temp = cumsum(y .- μ, dims = 1)
    return mean((maximum(temp, dims = 1) - minimum(temp, dims = 1)) / σ)
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


#=
count = 0
for δ in 1:5
    for x in 1:5
        for y in 1:5
            for z in 1:5
                count +=1
            end
        end
    end
end
count
4 * 5^4
2500 / 60
=#
