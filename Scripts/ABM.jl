using Random, Distributions
include("CoinTossXUtilities.jl")
clearconsole()
#---------------------------------------------------------------------------------------------------

#----- Structures -----#
struct Parameters
    Nᴸᵥ::Int64 # Number of fundamentalists for the low-frequency agent class
    Nᴸₜ::Int64 # Number of chartists for the low-frequency agent class
    Nᴴ::Int64 # Number of high-frequency agents
    λᴸmin::Float64 # Minimum inter-arrival rate parameter for low-frequency agents
    λᴸmax::Float64 # Maximum inter-arrival rate parameter for low-frequency agents
    λᴴmin::Float64 # Minimum inter-arrival rate parameter for high-frequency agents
    λᴴmax::Float64 # Maximum inter-arrival rate parameter for high-frequency agents
end
mutable struct Order
end
mutable struct LOBState
    sₜ::Int64 # Current spread
    ρₜ::Float64 # Current order imbalance
    mₜ::Float64 # Current mid-price
    bids::Dict{Int64, Order} # Stores all the active bids
    asks::Dict{Int64, Order} # Stores all the active asks
end
abstract type Agent end
mutable struct Chartist <: Agent
    dₜ::String # Decision to buy or sell
    xₜ₋₁::Float64 # Previous position
    xₜ::Float64 # Next position
    tradingGateway::Client # Object to submit orders
end
mutable struct Fundamentalist <: Agent
    dₜ::String # Decision to buy or sell
    fₜ::Float64 # Current perceived value
    xₜ₋₁::Float64 # Previous position
    xₜ::Float64 # Next position
    tradingGateway::Client # Object to submit orders
end
mutable struct HighFrequencyTrader

end
#---------------------------------------------------------------------------------------------------

#----- Initialize agents -----#
function InitializeAgents(parameters::Parameters)
    agents = Vector{Agent}()
    for i in 1:parameters.N
        # Chartists
        push!(agents, Chartist(...))
        # Fundamentalists
        push!(agents, Fundamentalist(...))
    end
    return chartists, fundamentalists
end
#---------------------------------------------------------------------------------------------------

#----- Update positions -----#
function UpdatePosition!(agent::HighFrequencyTrader, LOBState::LOBState)

end
function UpdatePosition!(agent::Fundamentalist, LOBState::LOBState)
    agent.dₜ = agent.fₜ < LOBState.mₜ ? "Sell" : "Buy"
end
function UpdatePosition!(agent::Chartist, LOBState::LOBState)

end
#---------------------------------------------------------------------------------------------------

#----- Update LOB state -----#
function UpdateLOBState!(LOB::LOBState, message::String)
    fields = split(message, ",")
    type = Symbol(fields[2]); side = Symbol(fields[3]); trader = Symbol(fields[4])
    executions = split(fields[5][2:(end - 1)], ",")
    for execution in executions
        executionFields = split(execution[2:(end - 1)])
        id = Int(executionFields[1]); price = Int(executionFields[2]); volume = Int(executionFields[3])
        if type == :New
            side == :Buy ? push!(LOB.bids, Order(id, price, volume, trader)) : push!(LOB.asks, Order(id, price, volume, trader))
        elseif type == :Cancelled
            side == :Buy ? delete!(LOB.bids, id) : push!(LOB.asks, id)
        elseif type == :Trade
            if side == :Buy
                LOB.asks[id] -= volume
                if LOB.asks[id] == 0
                    delete!(LOBState.asks, id)
                end
            else
                LOB.bids[id] -= volume
                if LOB.bids[id] == 0
                    delete!(LOBState.bids, id)
                end
            end
        end
    end
    bestBid = maximum(order -> order.price, values(LOBState.bids)); bestAsk = minimum(order -> order.price, values(LOBState.asks))
    totalBuyVolume = sum(order.volume for order in values(LOBState.bids)); totalSellVolume = sum(order.volume for order in values(LOBState.asks))
    LOBState.sₜ = abs(bestAsk - bestBid)
    LOBState.mₜ = (bestAsk + bestBid) / 2
    LOBState.ρₜ = (totalBuyVolume - totalSellVolume) / (totalBuyVolume + totalSellVolume)
    @spawnat Dosomethingwitheachagentclass
    @spawnat Dosomethingwitheachagentclass
end

#---------------------------------------------------------------------------------------------------

#----- Simulation -----#
function Simulate(horizon::Int64, parameters::Parameters; seed = 1)
    Random.seed!(seed)
    StartCoinTossX()
    client = Login(1, 1)
    try # This ensures that the client gets logged out whether an error occurs or not
        agents = InitializeAgents(parameters) # Initialize agents
        Juno.progress() do id # Progress bar
            t₀ = time()
            while (t - t₀) < horizon
                for agent in agents # Iterate through all agents
                    UpdatePosition!(agent, LOB)
                end
                t = time()
                @info "Trading" progress=((t - t₀) / horizon) _id=id # Update progress
            end
        end
    finally
        Logout(client)
    end
end
#---------------------------------------------------------------------------------------------------

#----- Supplementary functions -----#
function PowerLaw(xₘ, α) # Volumes
    return xₘ / (rand() ^ (1 / α))
end
#---------------------------------------------------------------------------------------------------
