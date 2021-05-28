#=
DataCleaning:
- Julia version: 1.5.3
- Authors: Ivan Jericevich, Patrick Chang, Tim Gebbie
- Function: Clean CoinTossX simulated data into L1LOB data for stylized fact analysis as well as visualisation of HFT time-series
- Structure:
    1. Supplementary functions
    2. Clean raw data into L1LOB format
    3. Plot simulation results
- Examples:
    taqData = CleanData("TAQ")
    VisualiseSimulation(taqData, "L1LOB"; format = "png")
=#
using CSV, DataFrames, Dates, Plots
clearconsole()
#---------------------------------------------------------------------------------------------------

#----- Supplementary functions -----#
mutable struct Order
    Price::Int64
    Volume::Int64
end
mutable struct Best
    Price::Int64
    Volume::Int64
    IDs::Vector{Int64}
end
function MidPrice(best::Best, contraBest::Best)
    return (isempty(best.IDs) || isempty(contraBest.IDs)) ? missing : (best.Price + contraBest.Price) / 2
end
function MicroPrice(best::Best, contraBest::Best)
    return (isempty(best.IDs) || isempty(contraBest.IDs)) ? missing : (best.Price * best.Volume + contraBest.Price * contraBest.Volume) / (best.Volume + contraBest.Volume)
end
function Spread(best::Best, contraBest::Best)
    return (isempty(best.IDs) || isempty(contraBest.IDs)) ? missing : abs(best.Price - contraBest.Price)
end
function OrderImbalance(bids::Dict{Int64, Order}, asks::Dict{Int64, Order})
    if isempty(bids) && isempty(asks)
        return missing
    elseif isempty(bids)
        return -1
    elseif isempty(asks)
        return 1
    else
        totalBuyVolume = sum(order.Volume for order in values(bids))
        totalSellVolume = sum(order.Volume for order in values(asks))
        return (totalBuyVolume - totalSellVolume) / (totalBuyVolume + totalSellVolume)
    end
end
#---------------------------------------------------------------------------------------------------

#----- Clean raw data into L1LOB format -----#
#=
Function:
    - Update LOB and best with LO
    - Full crossed orders are also added to the LOB and then aggressed with subsequent effective MOs
    - Append mid-price, micro-price and spread info
Arguments:
    - file = output file to which L1LOB will be printed
    - order = order to be processed
    - best = best bid (ask) if bid (ask) LO
    - contraBest = best ask (bid) if ask (bid) LO
    - lob = bid (ask) side of LOB if bid (ask) LO
    - side = ∈ {-1, 1}
    - allowCrossing = should crossed orders be handled or not
Output:
    - Order id of crossed order (if any)
    - All other variables are updated in-place
=#
function ProcessLimitOrder!(file::IOStream, order::DataFrameRow, best::Best, contraBest::Best, lob::Dict{Int64, Order}, side::Int64, allowCrossing::Bool)
    crossedOrderId = nothing # Initialize. This resets when a new LO occurs
    if isempty(lob) || isempty(best.IDs) # If the dictionary is empty, this order automatically becomes best
        best.Price = order.Price; best.Volume = order.Volume; best.IDs = [order.OrderId]
        midPrice = MidPrice(best, contraBest); microPrice = MicroPrice(best, contraBest); spread = Spread(best, contraBest)
        println(file, string(order.DateTime, ",", best.Price, ",", best.Volume, ",Limit,", side, ",", midPrice, ",", microPrice, ",", spread))
    else # Otherwise find the best
        if (side * order.Price) > (side * best.Price) # Change best if price of current order better than the best (side == 1 => order.Price > best.Price) (side == -1 => order.Price < best.Price)
            best.Price = order.Price; best.Volume = order.Volume; best.IDs = [order.OrderId] # New best is created. This includes crossed orders since they are handled later when the EMOs are processed
            if !isempty(contraBest.IDs) && (side * order.Price) >= (side * contraBest.Price) # Crossing order => limit order becomes effective market order (to avoid an error first check if the otherside isn't empty)
                if allowCrossing # Either disallow crossing or treat crossed order as effective market order at the next iteration
                    println(string("Order ", order.OrderId, " crossed the spread"))
                    crossedOrderId = order.OrderId # Update the L1LOB with the crossed order. The trades occuring hereafter will be deducted from this as well as from the best on the contra side
                    order.Type = :CrossedLimit # Change the order type
                else
                    error("Negative spread at order " * string(order.OrderId))
                end
            else # Only print the LO if it is not a CLO. The CLO is only printed if it is partially filled after aggressing against the contra side
                midPrice = MidPrice(best, contraBest); microPrice = MicroPrice(best, contraBest); spread = Spread(best, contraBest)
                println(file, string(order.DateTime, ",", best.Price, ",", best.Volume, ",Limit,", side, ",", midPrice, ",", microPrice, ",", spread))
            end
        elseif order.Price == best.Price # Add the new order's volume and orderid to the best if they have the same price
            best.Volume += order.Volume; push!(best.IDs, order.OrderId) # Best is ammended by adding volume to best and appending the order id
            midPrice = MidPrice(best, contraBest); microPrice = MicroPrice(best, contraBest); spread = Spread(best, contraBest)
            println(file, string(order.DateTime, ",", best.Price, ",", best.Volume, ",Limit,", side, ",", midPrice, ",", microPrice, ",", spread))
        end
    end
    push!(lob, order.OrderId => Order(order.Price, order.Volume)) # New order is always pushed to LOB dictionary only after best is processed
    return crossedOrderId
end
#=
Function:
    - Update LOB and best with MO
    - Append mid-price, micro-price and spread info
Arguments:
    - file = output file to which L1LOB will be printed
    - order = order to be processed
    - best = best bid (ask) if bid (ask) LO
    - contraBest = best ask (bid) if ask (bid) LO
    - lob = bid (ask) side of LOB if bid (ask) LO
    - side = ∈ {-1, 1}
Output:
    - All variables are updated in-place
=#
function ProcessMarketOrder!(file::IOStream, order::DataFrameRow, best::Best, contraBest::Best, lob::Dict{Int64, Order}, side::Int64)
    contraOrder = lob[order.OrderId] # Extract order on contra side
    if order.Volume == best.Volume # Trade filled best - remove from LOB, and update best
        delete!(lob, order.OrderId) # Remove the order from the LOB
        if !isempty(lob) # If the LOB is non empty find the best
            bestPrice = side * maximum(x -> side * x.Price, values(lob)) # Find the new best price (bid => side == 1 so find max price) (ask => side == -1 so find min price)
            indeces = [k for (k, v) in lob if v.Price == bestPrice] # Find the order ids of the best
            best.Price = bestPrice; best.Volume = sum(lob[i].Volume for i in indeces); best.IDs = indeces # Update the best
        else # If the LOB is empty remove best
            best.Price = 0; best.Volume = 0; best.IDs = Vector{Int64}()
        end
    else # Trade partially filled best
        if order.Volume == contraOrder.Volume # Trade filled contra order - remove order from LOB, remove order from best, and update best
            delete!(lob, order.OrderId)
            best.Volume -= order.Volume; best.IDs = setdiff(best.IDs, order.OrderId)
        else # Trade partially filled contra order - update LOB, update best
            lob[order.OrderId].Volume -= order.Volume
            best.Volume -= order.Volume
        end
    end
    midPrice = MidPrice(best, contraBest); microPrice = MicroPrice(best, contraBest); spread = Spread(best, contraBest)
    !isempty(best.IDs) ? println(file, string(order.DateTime, ",", best.Price, ",", best.Volume, ",Limit,", side, ",", midPrice, ",", microPrice, ",", spread)) : println(file, string(order.DateTime, ",missing,missing,Limit,", side, ",missing,missing,missing"))
end
#=
Function:
    - Update LOB and best with OC
    - Append mid-price, micro-price and spread info
Arguments:
    - file = output file to which L1LOB will be printed
    - order = order to be processed
    - best = best bid (ask) if bid (ask) LO
    - contraBest = best ask (bid) if ask (bid) LO
    - lob = bid (ask) side of LOB if bid (ask) LO
    - side = ∈ {-1, 1}
Output:
    - All variables are updated in-place
=#
function ProcessCancelOrder!(file::IOStream, order::DataFrameRow, best::Best, contraBest::Best, lob::Dict{Int64, Order}, side::Int64)
    delete!(lob, order.OrderId) # Remove the order from the LOB
    if order.OrderId in best.IDs # Cancel hit the best
        if !isempty(lob) # Orders still remain in the LOB - find and update best
            bestPrice = side * maximum(x -> side * x.Price, values(lob)) # Find the new best price (bid => side == 1 so find max price) (ask => side == -1 so find min price)
            indeces = [k for (k,v) in lob if v.Price == bestPrice] # Find the order ids of the best
            best.Price = bestPrice; best.Volume = sum(lob[i].Volume for i in indeces); best.IDs = indeces # Update the best
            midPrice = MidPrice(best, contraBest); microPrice = MicroPrice(best, contraBest); spread = Spread(best, contraBest)
            println(file, string(order.DateTime, ",", best.Price, ",", best.Volume, ",Cancelled,", side, ",", midPrice, ",", microPrice, ",", spread))
        else # The buy side LOB was emptied - update best
            best.Price = 0; best.Volume = 0; best.IDs = Vector{Int64}()
            println(file, string(order.DateTime, ",missing,missing,Cancelled,", side, ",missing,missing,missing"))
        end
    end # OC did not hit best
end
#=
Function:
    - Update LOB and best with effective MO. This requires us to aggress the order against both side of the LOB
    - Append mid-price, micro-price and spread info
Arguments:
    - file = output file to which L1LOB will be printed
    - order = order to be processed
    - nextOrder = the order type of the next order to check if a partially filled crossed LO needs to be printed
    - best = best bid (ask) if bid (ask) LO
    - contraBest = best ask (bid) if ask (bid) LO
    - lob = bid (ask) side of LOB if bid (ask) LO
    - side = ∈ {-1, 1}
Output:
    - The order id of the crossed order if it hasn't been fully handled. Otherwise nothing
    - All other variables are updated in-place
=#
function ProcessEffectiveMarketOrder!(file::IOStream, order::DataFrameRow, nextOrder::Symbol, best::NamedTuple, contraBest::NamedTuple, lob::Dict{Int64, Order}, crossedOrderId::Int64, side::Int64) # Remove the executed quantity from the LO and push the remaining volume to the LOB
    if order.Volume == best.Volume # Crossed LO (which is now at the best) was fully executed against contra side - remove from LOB, and update best. Note that all crossed orders will sit at the best
        delete!(lob, crossedOrderId) # Remove the order from the LOB
        if !isempty(lob) # If the LOB is non empty find the best
            bestPrice = side * maximum(x -> side * x.Price, values(lob)) # Find the new best price (bid => side == 1 so find max price) (ask => side == -1 so find min price)
            indeces = [k for (k,v) in lob if v.Price == bestPrice] # Find the order ids of the best
            best.Price = bestPrice; best.Volume = sum(lob[i].Volume for i in indeces); best.IDs = indeces # Update the best
        else # If the LOB is empty remove best
            best.Price = 0; best.Volume = 0; best.IDs = Vector{Int64}()
        end
        crossedOrderId = nothing # Crossed LO was filled so reset
    else # Trade partially filled best
        lob[crossedOrderId].Volume -= order.Volume # Update the volume of the crossed LO in the LOB
        best.Volume -= order.Volume
        if nextOrder != :Trade # If the crossed LO has been partially filled but aggressed all contra best orders (no more EMOs left)
            midPrice = MidPrice(best, contraBest); microPrice = MicroPrice(best, contraBest); spread = Spread(best, contraBest)
            println(file, string(order.DateTime, ",", best.Price, ",", best.Volume, ",CrossedLimit,", side, ",", midPrice, ",", microPrice, ",", spread)) # The best already updated so all thats needed is to print and set crossed order = nothing
            crossedOrderId = nothing # Reset the crossed order indicator as it has been fully handled (need not update crossed order indicator elsewhere)
        end
    end
    order.Type = :EffectiveMarket # Change the order type
    return crossedOrderId # Update the crossed order indicator
end
#=
Function:
    - Process all orders and clean raw TAQ data into L1LOB bloomberg format
Arguments:
    - orders = TAQ data
    - allowCrossing = should crossed orders be handled or not
Output:
    - TAQ data
    - Output L1LOB file written to csv
=#
function CleanData(taqFile::String; allowCrossing::Bool = false)
    orders = CSV.File(string("Data/", taqFile, ".csv"), types = Dict(:ClientOrderId => Int64, :DateTime => DateTime, :Price => Int64, :Volume => Int64, :Side => Symbol, :Type => Symbol, :TraderMnemonic => Int64), dateformat = "yyyy-mm-dd HH:MM:SS.s") |> DataFrame
    orders.Type[findall(x -> x == :New, orders.Type)] .= :Limit; orders.Type[findall(x -> x == 0, orders.Price)] .= :Market # Rename Types
    rename!(orders, [:ClientOrderId => :OrderId, :TraderMnemonic => :Trader])
    traders = CSV.File("Data/Trader.csv", drop = [:TraderId]) |> Tables.matrix |> vec
    orders.Trader = traders[orders.Trader] # Map tarder IDs to their label
    bids = Dict{Int64, Order}(); asks = Dict{Int64, Order}() # Both sides of the entire LOB are tracked with keys corresponding to orderIds
    bestBid = Best(0, 0, Vector{Int64}()); bestAsk = Best(0, 0, Vector{Int64}()) # Current best bid/ask is stored in a tuple (Price, vector of Volumes, vector of OrderIds) and tracked
    crossedOrderId = nothing # Stores the order that crossed if it hasn't been fully handled yet. If it has been processed then contains nothing
    imbalance = Vector{Union{Missing, Float64}}() # Calculate the order imbalance in the LOB
    open("Data/L1LOB.csv", "w") do file
        println(file, "DateTime,Price,Volume,Type,Side,MidPrice,MicroPrice,Spread") # Header
        Juno.progress() do id # Progress bar
            for i in 1:nrow(orders) # Iterate through all orders
                order = orders[i, :]
                #-- Limit Orders --#
                if order.Type == :Limit
                    if order.Side == :Buy # Buy limit order
                        crossedOrderId = ProcessLimitOrder!(file, order, bestBid, bestAsk, bids, 1, allowCrossing) # Add the order to the lob and update the best if necessary (crossed orders as well)
                    else # Sell limit order
                        crossedOrderId = ProcessLimitOrder!(file, order, bestAsk, bestBid, asks, -1, allowCrossing) # Add the order to the lob and update the best if necessary (crossed orders as well)
                    end
                #-- Market Orders --#
                elseif order.Type == :Trade # Market order always affects the best
                    if order.Side == :Sell # Trade was buyer-initiated (Sell MO)
                        println(file, string(order.DateTime, ",", order.Price, ",", order.Volume, ",Market,-1,missing,missing,missing")) # Sell trade is always printed (EMOs are just classified as MOs)
                        ProcessMarketOrder!(file, order, bestBid, bestAsk, bids, 1) # Sell trade affects bid side. Always aggress MO/EMO against contra side and update LOB and best. This is done before the processing of the crossed LO
                        if !isnothing(crossedOrderId) # The crossed order hasn't been fully handled yet. So after handling the above EMO aggress it against the crossed LO
                            crossedOrderId = ProcessEffectiveMarketOrder!(file, order, orders[i + 1, :Type], bestAsk, bestBid, asks, crossedOrderId, -1) # Aggress effective sell MO against ask side (originating crossed LO). Update crossed order at LOB and best. Print partially filled crossed LO if necessary
                        end
                    else # Trade was seller-initiated (Buy MO)
                        println(file, string(order.DateTime, ",", order.Price, ",", order.Volume, ",Market,1,missing,missing,missing")) # Buy trade is always printed (EMOs are just classified as MOs)
                        ProcessMarketOrder!(file, order, bestAsk, bestBid, asks, -1) # Buy trade affects ask side. Always aggress MO/EMO against contra side and update LOB and best. This is done before the processing of the crossed LO
                        if !isnothing(crossedOrderId) # The crossed order hasn't been fully handled yet. So after handling the above EMO aggress it against the crossed LO
                            crossedOrderId = ProcessEffectiveMarketOrder!(file, order, orders[i + 1, :Type], bestBid, bestAsk, bids, crossedOrderId, 1) # Aggress effective buy MO against bid side (originating crossed LO). Update crossed order at  LOB and best. Print partially filled crossed LO if necessary
                        end
                    end
                #-- Cancel Orders --#
                elseif order.Type == :Cancelled
                    if order.Side == :Buy # Cancel buy limit order
                        ProcessCancelOrder!(file, order, bestBid, bestAsk, bids, 1) # Aggress cancel order against buy side and update LOB and best
                    else # Cancel sell limit order
                        ProcessCancelOrder!(file, order, bestAsk, bestBid, asks, -1) # Aggress cancel order against sell side and update LOB and best
                    end
                end
                push!(imbalance, OrderImbalance(bids, asks)) # Calculate the volume imbalance after every iteration
                @info "Cleaning:" progress=(i / nrow(orders)) _id=id # Update progress
            end
        end
    end
    Juno.notification("Data cleaning complete"; kind = :Info, options = Dict(:dismissable => false))
    orders.DateTime = Dates.value.(orders.DateTime .- orders.DateTime[1]) ./ 1000 # Reformat times into relative times (milliseconds)
    orders.Imbalance = imbalance # Append order imbalance to data
    return orders
end
#---------------------------------------------------------------------------------------------------

#----- Plot simulation results -----#
function VisualiseSimulation(orders::DataFrame, l1lob::String; format = "pdf", duration = missing)
    # Cleaning
    filter!(x -> x.Type != :Market, orders)
    l1lob = CSV.File(string("Data/", l1lob, ".csv"), drop = [:Price, :Volume, :Side], missingstring = "missing") |> DataFrame |> x -> filter(y -> x.Type != :Market, x) # Filter out trades from L1LOB since their mid-prices are missing
    l1lob.DateTime = Dates.value.(l1lob.DateTime .- l1lob.DateTime[1]) ./ 1000
    if !ismissing(duration)
        filter!(x -> x.DateTime <= duration, l1lob); filter!(x -> x.DateTime <= duration, orders)
    end
    asks = filter(x -> x.Type == :Limit && x.Side == :Sell, orders); bids = filter(x -> x.Type == :Limit && x.Side == :Buy, orders)
    sells = filter(x -> x.Type == :Trade && x.Side == :Buy, orders); buys = filter(x -> x.Type == :Trade && x.Side == :Sell, orders) # scale * (asks.Volume / mean(asks.Volume))
    cancelAsks = filter(x -> x.Type == :Cancelled && x.Side == :Sell, orders); cancelBids = filter(x -> x.Type == :Cancelled && x.Side == :Buy, orders)
    # Bubble plot
    bubblePlot = plot(asks.DateTime, asks.Price, seriestype = :scatter, marker = (:red, stroke(:red), 0.7), label = "Ask (LO)", ylabel = "Price (ticks)", legend = :topleft, legendfontsize = 5, xrotation = 30)
    plot!(bubblePlot, bids.DateTime, bids.Price, seriestype = :scatter, marker = (:blue, stroke(:blue), 0.7), label = "Bid (LO)")
    plot!(bubblePlot, sells.DateTime, sells.Price, seriestype = :scatter, marker = (:red, stroke(:red), :utriangle, 0.7), label = "Sell (MO)")
    plot!(bubblePlot, buys.DateTime, buys.Price, seriestype = :scatter, marker = (:blue, stroke(:blue), :utriangle, 0.7), label = "Buy (MO)")
    plot!(bubblePlot, cancelAsks.DateTime, cancelAsks.Price, seriestype = :scatter, marker = (:red, stroke(:red), :xcross, 0.7), label = "Cancel Ask")
    plot!(bubblePlot, cancelBids.DateTime, cancelBids.Price, seriestype = :scatter, marker = (:blue, stroke(:blue), :xcross, 0.7), label = "Cancel Bid")
    # L1LOB features
    plot!(bubblePlot, l1lob.DateTime, l1lob.MidPrice, seriestype = :steppost, linecolor = :black, label = "Mid-price")
    plot!(bubblePlot, l1lob.DateTime, l1lob.MicroPrice, seriestype = :line, linecolor = :green, label = "Micro-price")
    # Spread and imbalance features
    volumeImbalance = plot(orders.DateTime, orders.Imbalance, seriestype = :line, linecolor = :purple, xlabel = "Time (s)", ylabel = "Order Imbalance", label = "OrderImbalance", legend = :topleft, legendfontsize = 5)
    plot!(twinx(), l1lob.DateTime, l1lob.Spread, seriestype = :steppost, linecolor = :orange, xlabel = "Time (s)", ylabel = "Spread", label = "Spread", legend = :topright, legendfontsize = 5)
    l = @layout([a; b{0.3h}])
    simulation = plot(bubblePlot, volumeImbalance, layout = l, link = :x, guidefontsize = 7)
    savefig(bubblePlot, "Figures/Simulation." * format)
    # Agents
    TFSells = filter(x -> x.Trader[1:2] == "TF", asks); TFBuys = filter(x -> x.Trader[1:2] == "TF", bids)
    TFAgents = plot(TFSells.DateTime, TFSells.Price, seriestype = :scatter, marker = (:red, stroke(:red), 0.7), label = "Sell (MO)", ylabel = "Price (ticks)", legend = :topleft, legendfontsize = 5, xrotation = 30)
    plot!(TFAgents, TFBuys.DateTime, TFBuys.Price, seriestype = :scatter, marker = (:blue, stroke(:blue), 0.7), label = "Buy (MO)")
    savefig(TFAgents, "Figures/TFAgents." * format)
    VISells = filter(x -> x.Trader[1:2] == "VI", sells); VIBuys = filter(x -> x.Trader[1:2] == "VI", buys)
    VIAgents = plot(VISells.DateTime, VISells.Price, seriestype = :scatter, marker = (:red, stroke(:red), :utriangle, 0.7), label = "Sell (MO)")
    plot!(VIAgents, VIBuys.DateTime, VIBuys.Price, seriestype = :scatter, marker = (:blue, stroke(:blue), :utriangle, 0.7), label = "Buy (MO)")
    savefig(VIAgents, "Figures/VIAgents." * format)
end
#---------------------------------------------------------------------------------------------------
