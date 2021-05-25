using Random, Distributions, Reactive, DataFrames
include("CoinTossXUtilities.jl")
clearconsole()
#---------------------------------------------------------------------------------------------------

#----- Structures -----#
struct Parameters
    Nᴸₜ::Int64 # Number of chartists for the low-frequency agent class
    Nᴸᵥ::Int64 # Number of fundamentalists for the low-frequency agent class
    Nᴴ::Int64 # Number of high-frequency agents
    λᴸ::Float64 # Inter-arrival rate parameter for low-frequency agents
    λᴴ::Float64 # Inter-arrival rate parameter for high-frequency agents
    λᴸmin::Float64 # Minimum inter-arrival rate parameter for low-frequency agents
    λᴸmax::Float64 # Maximum inter-arrival rate parameter for low-frequency agents
    λᴴmin::Float64 # Minimum inter-arrival rate parameter for high-frequency agents
    λᴴmax::Float64 # Maximum inter-arrival rate parameter for high-frequency agents
    δ::Float64     # Upper cut-off for LF agents decision rule
    T::Float64  # Simulation time
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
    p̄ₜ::Float64    # Agent's mid-price EWMA
    ActionTimes::Array{Float64,1} # Arrival process of when each agent makes a decision
    τ::Float64    # Time constant (used for EWMA computation) = Average inter-arrival of their decision times
end
mutable struct Fundamentalist <: Agent
    fₜ::Float64 # Current perceived value
    ActionTimes::Array{Float64,1} # Arrival process of when each agent makes a decision
end
mutable struct HighFrequency <: Agent
    ActionTimes::Array{Float64,1} # Arrival process of when each agent makes a decision
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

#----- Initialize agents and their parameters -----#
function InitializeAgents(parameters::Parameters)
    HFagents = Vector{Agent}()
    chartists = Vector{Agent}()
    fundamentalists = Vector{Agent}()
    for i in 1:parameters.Nᴴ
        push!(HFagents, HighFrequency(AgentTimes(parameters.λᴸ, parameters.λᴸmin, parameters.λᴸmax, parameters.T)))
    end
    for i in 1:parameters.Nᴸₜ
        ActionTimes = AgentTimes(parameters.λᴸ, parameters.λᴸmin, parameters.λᴸmax, parameters.T)
        push!(chartists, Chartist(100, ActionTimes, convert(Float64, mean(diff(ActionTimes)))))
    end
    for i in 1:parameters.Nᴸᵥ
        push!(fundamentalists, Fundamentalist(100, AgentTimes(parameters.λᴸ, parameters.λᴴmin, parameters.λᴴmax, parameters.T)))
    end
    return HFagents, chartists, fundamentalists
end
#---------------------------------------------------------------------------------------------------

#----- Agent rules -----#
function HighFrequencyAgentAction(LOB::LOBState, parameters::Parameters)
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
function ChartistAction(LOB::LOBState, chartists::Vector{Agent}, parameters::Parameters, agentNumber::Int64, orderNumber::Int64)
    # Update the agent's EWMA
    Δt = chartists[agentNumber].ActionTimes[orderNumber+1] - chartists[agentNumber].ActionTimes[orderNumber]    # Inter-arrival
    λ = 1 - exp(- Δt / chartists[agentNumber].τ)    # Decay for low-pass filter
    chartists[agentNumber].p̄ₜ += λ * (LOB.mₜ - chartists[agentNumber].p̄ₜ)
    # Check if agent performs action
    if abs(LOB.mₜ - chartists[agentNumber].p̄ₜ) <= LOB.sₜ
        return nothing
    end
    # Determine the lower volume bound based on agent decision rule
    if LOB.sₜ < abs(LOB.mₜ - chartists[agentNumber].p̄ₜ) <= (parameters.δ * LOB.mₜ)
        xₘ = 20
    end
    if abs(LOB.mₜ - chartists[agentNumber].p̄ₜ) > (parameters.δ * LOB.mₜ)
        xₘ = 50
    end
    decision = LOB.mₜ < chartists[agentNumber].p̄ₜ ? :Ask : :Bid
    α = decision == :Ask ? 1 - LOB.ρₜ : 1 + LOB.ρₜ
    volume = PowerLaw(xₘ, α)
    return decision, volume
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
# Initialize agents and parameters
parameters = Parameters(10, 10, 50, 20, 15, 10, 40, 10, 40, 0.05, 500)
# Initialize the agents
(HFagents, chartists, fundamentalists) = InitializeAgents(parameters)
# Initialize decision times
function CreateAgentDecisions(parameters::Parameters, HFagents::Vector{Agent}, chartists::Vector{Agent}, fundamentalists::Vector{Agent})
    # Create the dataframe
    DecisionTimes = DataFrame(Times = Float64[], OrderType = :Symbol, AgentType = :Symbol, AgentNumber = Int64[], AgentOrderNumber = Int64[])
    # Adding HF agents decsion and cancellation times to DecisionTimes
    for i in 1:parameters.Nᴴ
        ActionTimes = HFagents[i].ActionTimes[2:end]
        CancelTimes = filter(x -> x < parameters.T, ActionTimes .+ 300)
        ActionDF = DataFrame(Times = ActionTimes, OrderType = :LO, AgentType = :HF, AgentNumber = Int(i), AgentOrderNumber = collect(1:1:length(ActionTimes)))
        CancelDF = DataFrame(Times = CancelTimes, OrderType = :Cancel, AgentType = :HF, AgentNumber = Int(i), AgentOrderNumber = collect(1:1:length(CancelTimes)))
        DecisionTimes = [DecisionTimes; ActionDF]
        DecisionTimes = [DecisionTimes; CancelDF]
    end
    # Adding chartists decsion times to DecisionTimes
    for i in 1:parameters.Nᴸₜ
        ActionTimes = chartists[i].ActionTimes[2:end]
        ActionDF = DataFrame(Times = ActionTimes, OrderType = :MO, AgentType = :Chartist, AgentNumber = Int(i), AgentOrderNumber = collect(1:1:length(ActionTimes)))
        DecisionTimes = [DecisionTimes; ActionDF]
    end
    # Adding fundamentalists decsion times to DecisionTimes
    for i in 1:parameters.Nᴸᵥ
        ActionTimes = fundamentalists[i].ActionTimes[2:end]
        ActionDF = DataFrame(Times = ActionTimes, OrderType = :MO, AgentType = :Fundamentalist, AgentNumber = Int(i), AgentOrderNumber = collect(1:1:length(ActionTimes)))
        DecisionTimes = [DecisionTimes; ActionDF]
    end
    sort!(DecisionTimes, [:Times])
    return DecisionTimes
end
DecisionTimes = CreateAgentDecisions(parameters, HFagents, chartists, fundamentalists)

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
function rexp(θ, θmin, θmax)
    return rand(Distributions.truncated(Exponential(θ), θmin, θmax))
end
function AgentTimes(θ, θmin, θmax, T)
    t = Float64[]
    counter = 0.0
    push!(t, counter)
    while true
        τ = rexp(θ, θmin, θmax)
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
