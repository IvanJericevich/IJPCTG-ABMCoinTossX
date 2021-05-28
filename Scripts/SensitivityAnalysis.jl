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
