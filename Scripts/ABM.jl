using Random, Distributions, DataFrames, Dates, Sockets
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
    δ::Float64 # Upper cut-off for LF agents decision rule
    κ::Float64 # Scaling factor for order placement depth
    ν::Float64 # Scaling factor for power law order size
    m₀::Float64 # Initial mid-price
    σ::Float64 # Std dev in Normal for log-normal for fundamental value
    T::Millisecond # Simulation time
end
mutable struct LimitOrder
    price::Int64
    volume::Int64
    trader::Symbol
end
mutable struct MovingAverage
    p̄ₜ::Float64 # General EWMA of mid-price
    actionTimes::Array{Millisecond,1} # Arrival process of when all agents makes a decision
    τ::Float64 # Time constant (used for EWMA computation) = Average inter-arrival of their decision times
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
    p̄ₜ::Float64 # Agent's mid-price EWMA
    actionTimes::Array{Millisecond,1} # Arrival process of when each agent makes a decision
    τ::Float64 # Time constant (used for EWMA computation) = Average inter-arrival of their decision times
	orderNumber::Int64
end
mutable struct Fundamentalist <: Agent
    fₜ::Float64 # Current perceived value
    actionTimes::Array{Millisecond,1} # Arrival process of when each agent makes a decision
end
mutable struct HighFrequency <: Agent
    actionTimes::Array{Millisecond,1} # Arrival process of when each agent makes a decision
end
#---------------------------------------------------------------------------------------------------

#----- Initialize agents and their parameters -----#
function InitializeAgents(parameters::Parameters)
	chartists = Vector{Chartist}()
	HFagents = map(i -> HighFrequency(AgentTimes(parameters.λᴸ, parameters.λᴸmin, parameters.λᴸmax, parameters.T)), 1:parameters.Nᴴ)
	for i in 1:parameters.Nᴸₜ
		actionTimes = AgentTimes(parameters.λᴸ, parameters.λᴸmin, parameters.λᴸmax, parameters.T)
		push!(chartists, Chartist(parameters.m₀, actionTimes, convert(Float64, mean(diff(Dates.value.(actionTimes)))), 1))
	end
	fundamentalists = map(i -> Fundamentalist(parameters.m₀ * exp(rand(Normal(0, parameters.σ))), AgentTimes(parameters.λᴸ, parameters.λᴴmin, parameters.λᴴmax, parameters.T)), 1:parameters.Nᴸᵥ)
    return HFagents, chartists, fundamentalists
end
#---------------------------------------------------------------------------------------------------

#----- Agent rules -----#
function HighFrequencyAgentAction!(order::Order, LOB::LOBState, parameters::Parameters)
	if order.type == "Limit"
		θ = LOB.ρₜ/2 + .5 # Probability of placing an ask
	    order.side = rand() < θ ? "Sell" : "Buy"
	    if order.side == "Sell"
	        α = 1 - (LOB.ρₜ/parameters.ν) # Shape for power law
            η = floor(rand(Gamma(LOB.sₜ, exp(LOB.ρₜ / parameters.κ))))
	        order.price = LOB.bₜ + 1 + η
	        order.volume = round(Int, PowerLaw(10, α))
            order.displayVolume = order.volume
	    else
	        α = 1 + (LOB.ρₜ/parameters.ν)
            η = floor(rand(Gamma(LOB.sₜ, exp(-LOB.ρₜ / parameters.κ))))
	        order.price = LOB.aₜ - 1 - η
	        order.volume = round(Int, PowerLaw(10, α))
            order.displayVolume = order.volume
	    end
	else
		if haskey(LOB.bids, order.orderId)
			order.side = "Buy"
			order.price = LOB.bids[order.orderId].price
		elseif haskey(LOB.asks, order.orderId)
			order.side = "Sell"
			order.price = LOB.asks[order.orderId].price
		end
	end
end
function ChartistAction!(order::Order, LOB::LOBState, chartist::Chartist, parameters::Parameters)
    Δt = Dates.value(chartist.actionTimes[chartist.orderNumber + 1]) - Dates.value(chartist.actionTimes[chartist.orderNumber]) # Inter-arrival
    λ = 1 - exp(- Δt / chartist.τ) # Decay for low-pass filter
    chartist.p̄ₜ += λ * (LOB.mₜ - chartist.p̄ₜ)
    if !(abs(LOB.mₜ - chartist.p̄ₜ) <= LOB.sₜ) # Check if agent performs action
        # Determine the lower volume bound based on agent decision rule
        if LOB.sₜ < abs(LOB.mₜ - chartist.p̄ₜ) <= (parameters.δ * LOB.mₜ)
            xₘ = 20
        end
        if abs(LOB.mₜ - chartist.p̄ₜ) > (parameters.δ * LOB.mₜ)
            xₘ = 50
        end
        order.side = LOB.mₜ < chartist.p̄ₜ ? "Sell" : "Buy"
        α = order.side == "Sell" ? 1 - (LOB.ρₜ/parameters.ν) : 1 + (LOB.ρₜ/parameters.ν)
    end
	if (order.side == "Buy" && !isempty(LOB.asks)) || (order.side == "Sell" && !isempty(LOB.bids)) # Agent won't submit MO if no orders on contra side
		order.volume = round(Int, PowerLaw(xₘ, α))
	end
    # Update the agent's EWMA
    chartist.orderNumber += 1
end
function FundamentalistAction!(order::Order, LOB::LOBState, fundamentalists::Fundamentalist, parameters::Parameters)
    if !(abs(LOB.mₜ - fundamentalists.fₜ) <= LOB.sₜ) # Check if agent performs action
        # Determine the lower volume bound based on agent decision rule
        if LOB.sₜ < abs(LOB.mₜ - fundamentalists.fₜ) <= (parameters.δ * LOB.mₜ)
            xₘ = 20
        end
        if abs(LOB.mₜ - fundamentalists.fₜ) > (parameters.δ * LOB.mₜ)
            xₘ = 50
        end
        order.side = fundamentalists.fₜ < LOB.mₜ ? "Sell" : "Buy"
        α = order.side == "Sell" ? 1 - (LOB.ρₜ/parameters.ν) : 1 + (LOB.ρₜ/parameters.ν)
    end
	if (order.side == "Buy" && !isempty(LOB.asks)) || (order.side == "Sell" && !isempty(LOB.bids))
        order.volume = round(Int, PowerLaw(xₘ, α))
	end
end
#---------------------------------------------------------------------------------------------------

#----- Update Moving Average -----#
function UpdateMovingAverage!(LOB::LOBState, MA::MovingAverage, index::Int64)
    Δt = MA.actionTimes[index + 1] - MA.actionTimes[index] # Inter-arrival
    λ = 1 - exp(- Δt / MA.τ) # Decay for low-pass filter
    MA.p̄ₜ += λ * (LOB.mₜ - MA.p̄ₜ)
end
#---------------------------------------------------------------------------------------------------

#----- Update LOB state -----#
function UpdateLOBState!(LOB::LOBState, message)
	msg = split(message, "|")
	fields = split(msg[1], ",")
    type = Symbol(fields[1]); side = Symbol(fields[2]); trader = Symbol(fields[3])
    executions = msg[2:end]
	if executions[1] != ""
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
		totalBuyVolume = 0; totalSellVolume = 0
		if !isempty(LOB.bids) && !isempty(LOB.asks)
			LOB.bₜ = maximum(order -> order.price, values(LOB.bids)); LOB.aₜ = minimum(order -> order.price, values(LOB.asks))
			totalBuyVolume = sum(order.volume for order in values(LOB.bids)); totalSellVolume = sum(order.volume for order in values(LOB.asks))
		else
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
    t = Vector{Millisecond}()
    counter = Millisecond(0)
    push!(t, counter)
    while counter < T
        τ = Millisecond(round(Int, rexp(θ, θmin, θmax) * 1000))
        counter += τ
        push!(t, counter)
    end
    return t
end
#---------------------------------------------------------------------------------------------------

#----- Agent decision times -----#
function CreateAgentDecisions(parameters::Parameters, HFagents::Vector{HighFrequency}, chartists::Vector{Chartist}, fundamentalists::Vector{Fundamentalist}) # Initialize decision times
	decisionTimes = DataFrame(RelativeTime = Vector{Millisecond}(), Order = Vector{Order}(), AgentType = Vector{Symbol}(), AgentIndex = Vector{Int64}())
	idCounter = 0
    for i in 1:parameters.Nᴴ # Adding HF agents decsion and cancellation times to decisionTimes
        limitTimes = HFagents[i].actionTimes[2:end]
		limitOrders = map(x -> Order(orderId = idCounter + x, traderMnemonic = string("HF", i), type = "Limit"), 1:length(limitTimes))
		append!(decisionTimes, DataFrame(RelativeTime = limitTimes, Order = limitOrders, AgentType = :HF, AgentIndex = i)) # Limit orders
        cancelTimes = filter(x -> x < parameters.T, limitTimes .+ Millisecond(20 * 1000))
		cancelOrders = map(x -> Order(orderId = idCounter + x, traderMnemonic = string("HF", i), type = "Cancel"), 1:length(cancelTimes))
	    append!(decisionTimes, DataFrame(RelativeTime = cancelTimes, Order = cancelOrders, AgentType = :HF, AgentIndex = i))
		idCounter += length(limitTimes)
    end
    for i in 1:parameters.Nᴸₜ # Adding chartists decsion times to decisionTimes
        marketTimes = chartists[i].actionTimes[2:end]
		marketOrders = map(x -> Order(orderId = idCounter + x, traderMnemonic = string("TF", i), type = "Market"), 1:length(marketTimes))
        append!(decisionTimes, DataFrame(RelativeTime = marketTimes, Order = marketOrders, AgentType = :TF, AgentIndex = i))
		idCounter += length(marketTimes)
    end
    for i in 1:parameters.Nᴸᵥ # Adding fundamentalists decsion times to decisionTimes
		marketTimes = fundamentalists[i].actionTimes[2:end]
		marketOrders = map(x -> Order(orderId = idCounter + x, traderMnemonic = string("VI", i), type = "Market"), 1:length(marketTimes))
        append!(decisionTimes, DataFrame(RelativeTime = marketTimes, Order = marketOrders, AgentType = :VI, AgentIndex = i))
		idCounter += length(marketTimes)
    end
    sort!(decisionTimes, :RelativeTime)
    return decisionTimes
end
#---------------------------------------------------------------------------------------------------

#----- Simulation -----#
function InjectSimulation(parameters; startJVM = false, seed = 1)
	Random.seed!(seed)
    (HFagents, chartists, fundamentalists) = InitializeAgents(parameters) # Initialize the agents
    data = CreateAgentDecisions(parameters, HFagents, chartists, fundamentalists) # Initialize decision times
    # Initialize LOB state and MA
    LOB = LOBState(100, 0, 1000, 950, 1050, Dict{Int64, LimitOrder}(), Dict{Int64, LimitOrder}())
    MA = MovingAverage(parameters.m₀, [Millisecond(0); data.RelativeTime], mean(diff(Dates.value.([Millisecond(0); data.RelativeTime]))))
    times = Vector{Millisecond}(); midprice = Vector{Float64}() # Setup storage for time series of mid-price
    startJVM ? StartJVM() : nothing # Start JVM for first run, thereafter don't start it
    gateway = Login(1, 1)
    receiver = UDPSocket()
    connected = bind(receiver, ip"127.0.0.1", 1234)
    if connected
        println("Market data listener connected")
    end
    try # This ensures that the client gets logged out whether an error occurs or not
        data.Time = data.RelativeTime #.+ Time(now())
        Juno.progress() do id # Progress bar
            for i in 1:nrow(data)
				order = data[i, :Order]
                messagesent = false
				if data[i, :AgentType] == :HF # High-frequency
					HighFrequencyAgentAction!(order, LOB, parameters)
					# data[i, :Time] <= Time(now()) ? println(string("Timeout: ", Time(now()) - data[i, :Time])) : sleep(data[i, :Time] - Time(now()))
					if order.type == "Limit" # Limit order
	                    SubmitOrder(gateway, order)
                        messagesent = true
	                elseif order.type == "Cancel" && order.price != 0 # Order cancel
	                    CancelOrder(gateway, order)
                        messagesent = true
	                end
				elseif data[i, :AgentType] == :TF # Chartist
					ChartistAction!(order, LOB, chartists[data[i, :AgentIndex]], parameters)
					# data[i, :Time] <= Time(now()) ? println(string("Timeout: ", Time(now()) - data[i, :Time])) : sleep(data[i, :Time] - Time(now()))
					if order.volume != 0
						SubmitOrder(gateway, order)
                        messagesent = true
					end
				else # Fundamentalist
					FundamentalistAction!(order, LOB, fundamentalists[data[i, :AgentIndex]], parameters)
					# data[i, :Time] <= Time(now()) ? println(string("Timeout: ", Time(now()) - data[i, :Time])) : sleep(data[i, :Time] - Time(now()))
					if order.volume != 0
						SubmitOrder(gateway, order)
                        messagesent = true
					end
				end
                if messagesent == true
                    message = String(recv(receiver))
                    println(message)
                    UpdateLOBState!(LOB, message)
                    println(LOB.sₜ)
                    println(LOB.ρₜ)
                    push!(times, data[i, :Time])
                    push!(midprice, LOB.mₜ)
                end
                @info "Trading" progress=(data.RelativeTime[i] / data.RelativeTime[end]) _id=id # Update progress
            end
        end
    finally
        close(receiver)
        Logout(gateway)
    end
    return times, midprice
end
#---------------------------------------------------------------------------------------------------

#----- Example -----#
StartCoinTossX(build = false)
param = Parameters(10, 25, 15, 20, 20, 10, 40, 10, 40, 0.05, 1.5, 2, 1000, 0.5, Millisecond(3600 * 1000))
x = InjectSimulation(param; startJVM = true) # First run, must start JVM
x = InjectSimulation(param) # Next run must be on different security, don't need to start JVM again

using Plots
plot(Dates.value.(x[1]) ./ 1000, x[2])

StopCoinTossX()
exit()
#---------------------------------------------------------------------------------------------------
