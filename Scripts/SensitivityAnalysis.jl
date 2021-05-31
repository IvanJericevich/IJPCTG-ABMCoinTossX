using HypothesisTests, Statistics, StatsBase, GARCH

mutable struct Moments # Moments of log-returns
    μ::Float64 # Mean
    σ::Float64 # Standard deviation
    κ::Float64 # Kurtosis
    α::Float64 # Power-law parameter calibrated to ACF
    hurst::Float64 # Hurst exponent: hurst < 0.5 => mean reverting; hurst == 0.5 => random walk; hurst > 0.5 => momentum
    ks::Float64 # Kolmogorov-Smirnov statistic - quantifies a distance between the empirical distribution function of the sample and the cumulative distribution function of the reference distribution
    gph::Float64 # GPH estimator representing long-range dependence
    adf::Float64 # ADF statistic representing random walk property of returns
    garch::Float64 # GARCH paramaters representing short-range dependence
    hill::Float64 # Hill estimator
end

mean(returns)
std(returns)
kurtosis(returns)
ExactOneSampleKSTest(returns, Normal())
ADFTest(returns, deterministic, lag)

surfaces = Array{Moments, 2}()

δgrid
κgrid
νgrid
σgrid

δmoments = Dict{Float64, Array{Float64, 3}}() # Each key is a single value of δ
for δ in δgrid
    combinations = zeros(Float64, 2, 2, 2)
    for x in 1:length(κgrid)
        for y in 1:length(νgrid)
            for z in 1:length(σgrid)
                parameters = Parameters()
                time, midPrice = InjectSimulation()
                combinations[x, y, z] = ComputeMoments()
            end
        end
    end
    δmoments[δ] = combinations
end
5*4

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
