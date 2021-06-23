#=
Calibration:
- Julia version: 1.5.3
- Authors: Ivan Jericevich, Patrick Chang, Tim Gebbie
- Function: Calibrate CoinTossX ABM simulated mid-prices to JSE mid-price data
- Structure:
    1. Moving block bootstrap to estimate covariance matrix of empirical moments on JSE mid-price time-series
    2. Objective function to be minimized
    3. Calibrate with NMTA optimization
- Examples:
=#
include("CoinTossXUtilities.jl"); include("ABMVolatilityAuctionProxy.jl"); include("Moments.jl"); ; include("NMTA.jl")
#---------------------------------------------------------------------------------------------------

#----- Moving block bootstrap to estimate covariance matrix of empirical moments on JSE mid-price time-series -----#
function MovingBlockBootstrap()
    return W
end
#---------------------------------------------------------------------------------------------------

#----- Objective function to be minimized -----#
function WeightedSumofSquaredErrors(parameters::Parameters, W::Array{Float64, 2}, empiricalMoments::Moments, gateway::TradingGateway)
    simulation = InjectSimulation(gateway, parameters)
    simulatedMoments = Moments(simulation)
    error = simulatedMoments - empiricalMoments
    return transpose(error) * W * error
end
#---------------------------------------------------------------------------------------------------

#----- Calibrate with NMTA optimization -----#

#---------------------------------------------------------------------------------------------------
