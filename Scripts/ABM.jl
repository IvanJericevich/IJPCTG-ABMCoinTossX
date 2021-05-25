using Random, Distributions, Reactive
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
mutable struct LimitOrder
    price::Int64
    volume::Int64
    trader::Symbol
end
mutable struct LOBState
    sₜ::Int64 # Current spread
    ρₜ::Float64 # Current order imbalance
    mₜ::Float64 # Current mid-price
    bₜ::Int64 # Best bid
    aₜ::Int64 # Best ask
    bids::Dict{Int64, LimitOrder} # Stores all the active bids
    asks::Dict{Int64, LimitOrder} # Stores all the active asks
end
abstract type Agent end
mutable struct Chartist <: Agent
    dₜ::String # Decision to buy or sell
    tradingGateway::Client # Object to submit orders
end
mutable struct Fundamentalist <: Agent
    dₜ::String # Decision to buy or sell
    fₜ::Float64 # Current perceived value
    tradingGateway::Client # Object to submit orders
end
#---------------------------------------------------------------------------------------------------

#----- Market Data Listener -----#
function Listen(LOB::LOBState)
    receiver = UDPSocket()
    connected = bind(receiver, ip"127.0.0.1", 1234)
    if connected
		try
			while true
	            message = String(recv(receiver))
				UpdateLOBState!(LOB, message)
	        end
		finally
			close(receiver)
		end
    end
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

#----- Agent rules -----#
function HighFrequencyAgentAction(LOB::LOBState, agentNumber::Int64)
    θ = LOB.ρₜ/2 + .5   # Probability of placing an ask
    decision = rand() < θ ? :Ask : :Bid
    if decision == :Ask     # Sell order
        α = 1 - LOB.ρₜ  # Shape for power law
        λₜ = LOB.sₜ / exp(- LOB.ρₜ / 2)     # Placement depth parameter
        η = floor(- λₜ * log(rand()))
        limitPrice = LOB.bₜ + 1 + η
        volume = PowerLaw(10, α)
    else                    # Buy order
        α = 1 + LOB.ρₜ  # Shape for power law
        λₜ = LOB.sₜ / exp(LOB.ρₜ / 2)       # Placement depth parameter
        η = floor(- λₜ * log(rand()))
        limitPrice = LOB.aₜ - 1 - η
        volume = PowerLaw(10, α)
    end
    return decision, limitPrice, volume
end
#---------------------------------------------------------------------------------------------------

#----- Update positions -----#
function UpdatePosition!(agent::HighFrequencyTrader, LOB::LOBState)
end
function UpdatePosition!(agent::Fundamentalist, LOB::LOBState)
    agent.dₜ = agent.fₜ < LOB.mₜ ? "Sell" : "Buy"
end
function UpdatePosition!(agent::Chartist, LOB::LOBState)
end
#---------------------------------------------------------------------------------------------------

#----- Update LOB state -----#
function UpdateLOBState!(LOB::LOBState, message)
	msg = split(message, "|")
	fields = split(msg[1], ",")
    type = Symbol(fields[1]); side = Symbol(fields[2]); trader = Symbol(fields[3])
    executions = msg[2:end]
    for execution in executions
        executionFields = split(execution, ",")
        id = parse(Int, executionFields[1]); price = parse(Int, executionFields[2]); volume = parse(Int, executionFields[3])
        if type == :New
            side == :Buy ? push!(LOB.bids, id => LimitOrder(price, volume, trader)) : push!(LOB.asks, id => LimitOrder(price, volume, trader))
        elseif type == :Cancelled
            side == :Buy ? delete!(LOB.bids, id) : delete!(LOB.asks, id)
        elseif type == :Trade
            if side == :Buy
                LOB.asks[id].volume -= volume
                if LOB.asks[id].volume == 0
                    delete!(LOB.asks, id)
                end
            else
                LOB.bids[id].volume -= volume
                if LOB.bids[id].volume == 0
                    delete!(LOB.bids, id)
                end
            end
        end
    end
	totalBuyVolume = totalSellVolume = 0
	if !isempty(LOB.bids) && !isempty(LOB.asks) # If both sides are non-empty the update the spread and midprice
		LOB.bₜ = maximum(order -> order.price, values(LOB.bids)); LOB.aₜ = minimum(order -> order.price, values(LOB.asks))
		totalBuyVolume = sum(order.volume for order in values(LOB.bids)); totalSellVolume = sum(order.volume for order in values(LOB.asks))
	else # Otherwise only update the best and keep spread and mid-price as before
		if !isempty(LOB.bids)
			LOB.bₜ = maximum(order -> order.price, values(LOB.bids))
			totalBuyVolume = sum(order.volume for order in values(LOB.bids))
		end
		if !isempty(LOB.asks)
			LOB.aₜ = minimum(order -> order.price, values(LOB.asks))
			totalSellVolume = sum(order.volume for order in values(LOB.asks))
		end
	end
	LOB.sₜ = abs(LOB.aₜ - LOB.bₜ)
	LOB.mₜ = (LOB.aₜ + LOB.bₜ) / 2
    LOB.ρₜ = (totalBuyVolume == 0) && (totalSellVolume == 0) ? 0.0 : (totalBuyVolume - totalSellVolume) / (totalBuyVolume + totalSellVolume)
end
#---------------------------------------------------------------------------------------------------

#----- Preset agent decision times -----#
DecisionTimes = DataFrame(Times = Float64[], OrderType = :Symbol, AgentType = :Symbol, AgentNumber = Int64[])
T = 500
N𝐟 = 10
# Adding HF agents decsion and cancellation times to DecisionTimes
for i in 1:N𝐟
    ActionTimes = AgentTimes(15, T)
    CancelTimes = filter(x -> x<T, ActionTimes .+ 300)
    ActionDF = DataFrame(Times = ActionTimes, OrderType = :LO, AgentType = :HF, AgentNumber = Int(i))
    CancelDF = DataFrame(Times = CancelTimes, OrderType = :Cancel, AgentType = :HF, AgentNumber = Int(i))
    DecisionTimes = [DecisionTimes; ActionDF]
    DecisionTimes = [DecisionTimes; CancelDF]
end
sort!(DecisionTimes, [:Times])
#---------------------------------------------------------------------------------------------------

#----- Simulation -----#
function Simulate(horizon::Int64, parameters::Parameters; seed = 1)
    Random.seed!(seed)
    StartCoinTossX()
    client = Login(1, 1)
	LOB = LOBState(0, 0, 0, 0, 0, Dict{Int64, LimitOrder}(), Dict{Int64, LimitOrder}())
	@async Listen(LOB)
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
function rexp(mean)
    return -mean * log(rand())
end
function AgentTimes(mean, T)
    t = Float64[]
    counter = 0.0
    while true
        τ = rexp(mean)
        counter += τ
        if counter > T
            break
        end
        push!(t, counter)
    end
    return t
end
#---------------------------------------------------------------------------------------------------

#----- Example -----#
# StartCoinTossX()
# client = Login(1, 1)
# LOB = LOBState(0, 0, 0, 0, 0, Dict{Int64, LimitOrder}(), Dict{Int64, LimitOrder}())
# z = @async Listen(LOB)
# SubmitOrder(client, Order(1, "John", "Buy", "Limit", 100, 50))
# Logout(client)
# StopCoinTossX()
#---------------------------------------------------------------------------------------------------
