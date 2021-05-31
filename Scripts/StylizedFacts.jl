#=
StylizedFacts:
- Julia version: 1.5.3
- Authors: Ivan Jericevich, Patrick Chang, Dieter Hendricks, Tim Gebbie
- Function: Plot the stylized facts of HFT data for different time resolutions
- Structure:
    1. Log return sample distributions for different time resolutions
    2. Log-return and absolute log-return autocorrelation
    3. Trade sign autocorrealtion
    4. Trade inter-arrival time distribution
    5. Extreme log-return percentile distribution for different time resolutions
    6. Price Impact
- Examples
    LogReturnDistribution(:TickbyTick; lobFile = "Model2L1LOB", format = "png")
    LogReturnAutocorrelation(50, lobFile = "Model2L1LOB", format = "png")
    TradeSignAutocorrelation(20, lobFile = "Model2L1LOB", format = "png")
    InterArrivalTimeDistribution("Model2L1LOB"; format = "png")
    ExtremeLogReturnPercentileDistribution(:TickbyTick, :Upper; lobFile = "Model2L1LOB", format = "png")
    VPIN(data, 50, 10, Millisecond(60 * 30 * 1000))
- TODO:
- TODO: Insert plot annotations for the values of α when fitting power laws and excess kurtosis
- TODO: Change font sizes
=#
using Distributions, CSV, Plots, DataFrames, StatsPlots, Dates, StatsBase, LaTeXStrings
clearconsole()
#---------------------------------------------------------------------------------------------------

#----- Extract OHLCV data -----#
function OHLCV(lobFile::DataFrame, resolution)
    l1lob = CSV.File(string("Data/", lobFile, ".csv"), drop = [:Price, :Side, :Spread], missingstring = "missing") |> DataFrame
    barTimes = l1lob.DateTime[1]:resolution:l1lob.DateTime[end]
    open(string("Data/OHLCV.csv"), "w") do file
        println(file, "DateTime,MidOpen,MidHigh,MidLow,MidClose,MicroOpen,MicroHigh,MicroLow,MicroClose,Volume,VWAP")
        for t in 1:(length(barTimes) - 1)
            startIndex = searchsortedfirst(l1lob.DateTime, barTimes[t])
            endIndex = searchsortedlast(l1lob.DateTime, barTimes[t + 1])
            if !(startIndex >= endIndex)
                bar = l1lob[startIndex:endIndex, :]
                tradesBar = filter(x -> x.Type == :Market, bar)
                midPriceOHLCV = string(bar.MidPrice[1], ",", maximum(skipmissing(bar.MidPrice)), ",", minimum(skipmissing(bar.MidPrice)), ",", bar.MidPrice[end])
                microPriceOHLCV = string(bar.MicroPrice[1], ",", maximum(skipmissing(bar.MicroPrice)), ",", minimum(skipmissing(bar.MicroPrice)), ",", bar.MicroPrice[end])
                vwap = !isempty(tradesBar) ? sum(tradesBar.Volume .* tradesBar.Price) / sum(tradesBar.Volume) : missing
                println(file, string(barTimes[t], ",", midPriceOHLCV, ",", microPriceOHLCV, ",", sum(bar.Volume), ",", vwap))
            end
        end
    end
end
#---------------------------------------------------------------------------------------------------

#----- Log return sample distributions for different time resolutions -----#
function LogReturnDistribution(resolution::Symbol; lobFile::String, cummulative::Bool = false, format::String = "pdf")
    # Obtain log-returns of price series
    if resolution == :TickbyTick
        data = CSV.File(string("Data/", lobFile, ".csv"), missingstring = "missing") |> DataFrame
        logReturns = diff(log.(filter(x -> !ismissing(x), data[:, :MidPrice])))
    else
        data = CSV.File(string("Data/MicroPrice ", resolution, " Bars.csv"), missingstring = "missing") |> DataFrame
        logReturns = diff(log.(filter(x -> !ismissing(x), data[:, :Close])))
    end
    # Fit theoretical distributions
    theoreticalDistribution = fit(Normal, logReturns)
    if cummulative
        # Plot empirical distribution
        empiricalDistribution = bar(logReturns, func = cdf, fillcolor = :blue, linecolor = :blue, xlabel = "Log returns", ylabel = "Density", label = "Empirical", legendtitle = "Distribution")
        # Plot fitted theoretical distributions
        plot!(empiricalDistribution, cdf.(theoreticalDistribution, logReturns), linecolor = :black, label = "Fitted Normal")
        savefig(empiricalDistribution, string("Figures/Log-ReturnCummulativeDistribution", resolution, ".", format))
    else
        # Plot empirical distribution
        empiricalDistribution = histogram(logReturns, normalize = :pdf, fillcolor = :blue, linecolor = :blue, xlabel = "Log returns", ylabel = "Density", label = "Empirical", legendtitle = "Empirical Distribution", legend = :topleft)
        # Plot fitted theoretical distributions
        # plot!(empiricalDistribution, [theoreticalDistribution1, theoreticalDistribution2], linecolor = [:green :purple], label = ["Fitted Normal" "Fitted Inverse Gaussian"])
        plot!(empiricalDistribution, theoreticalDistribution, linecolor = :black, label = "Fitted Normal")
        qqplot!(empiricalDistribution, Normal, logReturns, xlabel = "Normal theoretical quantiles", ylabel = "Sample quantiles", linecolor = :black, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, marker = (:blue, stroke(:blue), 3), legend = false, inset = (1, bbox(0.65, 0.1, 0.33, 0.33, :top)), subplot = 2)
        savefig(empiricalDistribution, string("Figures/Log-ReturnDistribution", resolution, ".", format))
    end
end
#---------------------------------------------------------------------------------------------------

#----- Log-return and absolute log-return autocorrelation -----#
function LogReturnAutocorrelation(lag::Int64; lobFile::String, format::String = "pdf")
    # Obtain log-returns
    data = CSV.File(string("Data/", lobFile, ".csv"), missingstring = "missing") |> DataFrame
    logReturns = diff(log.(filter(x -> !ismissing(x), data.MicroPrice)))
    # Calculate autocorrelations
    autoCorr = autocor(logReturns, 1:lag; demean = false)
    absAutoCorr = autocor(abs.(logReturns), 1:lag; demean = false)
    # Plot
    autoCorrPlot = plot(autoCorr, seriestype = :sticks, linecolor = :black, legend = false, xlabel = "Lag", ylabel = "Autocorrelation")
    plot!(autoCorrPlot, [1.96 / sqrt(length(logReturns)), -1.96 / sqrt(length(logReturns))], seriestype = :hline, line = (:dash, :black, 1))
    plot!(autoCorrPlot, absAutoCorr, seriestype = :sticks, linecolor = :black, legend = false, xlabel = "Lag", ylabel = "Autocorrelation", guidefontsize = 7, inset = (1, bbox(0.62, 0.55, 0.33, 0.33, :top)), subplot = 2)
    savefig(autoCorrPlot, "Figures/Log-ReturnAutocorrelation." * format)
end
#---------------------------------------------------------------------------------------------------

#----- Trade sign autocorrealtion -----#
function TradeSignAutocorrelation(lag::Int64; lobFile::String, format::String = "pdf")
    # Extract trade signs
    data = CSV.File(string("Data/", lobFile, ".csv"), missingstring = "missing") |> DataFrame
    tradeSignsTemp = filter(x -> x.Type == "Trade" && x.Price == 0, data).Side
    tradeSigns = map(x -> x == "Buy" ? 1 : -1, tradeSignsTemp)
    # Calculate autocorrelations
    autoCorr = autocor(tradeSigns, 1:lag; demean = false)
    # Plot
    autoCorrPlot = plot(autoCorr, seriestype = :sticks, linecolor = :black, legend = false, xlabel = "Lag", ylabel = "Autocorrelation")
    plot!(autoCorrPlot, [quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(tradeSigns)), quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(tradeSigns))], seriestype = :hline, line = (:dash, :black, 1))
    plot!(autoCorrPlot, autoCorr, xscale = :log10, inset = (1, bbox(0.58, 0.0, 0.4, 0.4)), subplot = 2, legend = false, xlabel = "Lag", guidefontsize = 7, ylabel = "Autocorrelation", linecolor = :black) #  ", L"(\log_{10})
    savefig(autoCorrPlot, "Figures/Trade-SignAutocorrelation." * format)
end
#---------------------------------------------------------------------------------------------------

#----- Trade inter-arrival time distribution -----#
function InterArrivalTimeDistribution(lobFile::String; format::String = "pdf")
    # Extract inter-arrival times
    data = CSV.File(string("Data/", lobFile, ".csv"), missingstring = "missing") |> DataFrame
    interArrivals = filter(x -> x.Type == :MO, data) |> y -> diff(y.DateTime) |> z -> Dates.value.(round.(z, Dates.Millisecond))
    # Estimate power-law distribution parameters
    xₘᵢₙ = minimum(interArrivals)
    α = 1 + length(interArrivals) / sum(log.(interArrivals ./ xₘᵢₙ))
    # Extract theoretical quantiles
    theoreticalQuantiles = map(i -> (1 - (i / length(interArrivals))) ^ (-1 / (α - 1)) * xₘᵢₙ, 1:length(interArrivals))
    theoreticalDistribution1 = fit(Exponential, interArrivals)
    theoreticalDistribution2 = fit(Weibull, interArrivals)
    theoreticalDistribution3 = fit(LogNormal, interArrivals)
    # Plot
    logDistribution = histogram(interArrivals, normalize = :pdf, linecolor = :blue, fillcolor = :blue, yscale = :log10, xlabel = "Inter-arrival time (milliseconds)", ylabel = "Log Density", annotations = (3, y[3], Plots.text(string("α=" α), :left)), legend = false)
    plot!(logDistribution, [theoreticalDistribution1, theoreticalDistribution2, theoreticalDistribution3], linecolor = [:green :purple], label = ["Fitted Exponential" "Fitted Weibull" "Fitted Log-Normal"])
    plot!(logDistribution, [theoreticalQuantiles theoreticalQuantiles], [interArrivals theoreticalQuantiles], seriestype = [:scatter :line], inset = (1, bbox(0.6, 0.03, 0.34, 0.34, :top)), subplot = 2, legend = :none, xlabel = "Power-Law Theoretical Quantiles", ylabel = "Sample Quantiles", linecolor = :black, markercolor = :blue, markerstrokecolor = :blue, scale = :log10)
    savefig(logDistribution, "TradeLog-Inter-ArrivalDistribution." * format)
end
#---------------------------------------------------------------------------------------------------

#----- Extreme log-return percentile distribution for different time resolutions -----#
function ExtremeLogReturnPercentileDistribution(resolution::Symbol, side::Symbol; lobFile::String, format::String = "pdf")
    # Obtain log-returns of price series
    if resolution == :TickbyTick
        data = CSV.File(string("Data/", lobFile, ".csv"), missingstring = "missing") |> DataFrame
        logReturns = diff(log.(filter(x -> !ismissing(x), data[:, :MidPrice])))
    else
        data = CSV.File(string("MicroPrice ", resolution, " Bars.csv"), missingstring = "missing") |> DataFrame
        logReturns = diff(log.(filter(x -> !ismissing(x), data[:, :Close])))
    end
    sort!(logReturns)
    # Extract extreme empirical quantiles
    observations = side == :Upper ? logReturns[findall(x -> x >= quantile(logReturns, 0.95), logReturns)] : -logReturns[findall(x -> x <= quantile(logReturns, 0.05), logReturns)]
    # Estimate power-law distribution parameters
    xₘᵢₙ = minimum(observations)
    α = 1 + length(observations) / sum(log.(observations ./ xₘᵢₙ))
    # Extract theoretical quantiles
    theoreticalQuantiles = map(i -> (1 - (i / length(observations))) ^ (-1 / (α - 1)) * xₘᵢₙ, 1:length(observations))
    # Plot
    extremePercentileDistributionPlot = plot(observations, normalize = :pdf, linecolor = :blue, xlabel = string("Log return Extreme", side, " percentiles"), ylabel = "Density", legend = false)#, annotations = (3, y[3], Plots.text(string("α=" α), :left))))
    plot!(extremePercentileDistributionPlot, [theoreticalQuantiles theoreticalQuantiles], [observations theoreticalQuantiles], seriestype = [:scatter :line], inset = (1, bbox(0.2, 0.03, 0.34, 0.34, :top)), subplot = 2, guidefontsize = 7, legend = :none, xlabel = "Power-Law Theoretical Quantiles", ylabel = "Sample Quantiles", linecolor = :black, markercolor = :blue, markerstrokecolor = :blue, scale = :log10)
    savefig(extremePercentileDistributionPlot, string("Figures/Extreme", side, "Log-ReturnPercentilesDistribution", resolution,".", format))
end
#---------------------------------------------------------------------------------------------------

#----- VPIN -----#
function VPIN(data, V, n, resolution)
    filter!(x -> x.Type == :WalkingMO, data)
    timeIntervals = collect(data.DateTime[1]:resolution:data.DateTime[end])
    VPIN = fill(0.0, length(timeIntervals) - 1)
    expandedData = Vector{Tuple{Millisecond, Symbol}}()
    for trade in eachrow(data)
        for v in 1:trade.Volume
            push!(expandedData, (trade.DateTime, trade.Side))
        end
    end
    for t in 2:length(timeIntervals)
        unitTrades = last.(expandedData[findall(x -> timeIntervals[t - 1] <= first(x) < timeIntervals[t], expandedData)])
        Vₜᴮ = Vector{Int64}(); Vₜᴬ = Vector{Int64}()
        L = 0; τ = 0; I = length(unitTrades) # Initialisation
        while true
            τ += 1
            if I >= (τ * V)
                push!(Vₜᴮ, sum(unitTrades .== :Buy))
                push!(Vₜᴬ, sum(unitTrades .== :Sell))
            else
                L = τ - 1
                break
            end
        end
        if L >= n
            VPIN[t - 1] = sum(abs.(Vₜᴬ[(L - n + 1):L] - Vₜᴮ[(L - n + 1):L])) / (n * V)
        else
            VPIN[t - 1] = NaN
        end
    end
    return VPIN
end
#---------------------------------------------------------------------------------------------------

#----- Extract price-impact data -----#
function PriceImpact(file::String)
    data = CSV.File(string("Data/", file, ".csv"), drop = [:MicroPrice, :Spread, :Imbalance], types = Dict(:Type => Symbol, :Price => Int64, :Volume => Int64, :MidPrice => Float64), missingstring = "missing") |> DataFrame
    tradeIndeces = findall(x -> x == :MO, data.Type); totalTradeCount = length(tradeIndeces)
    dayVolume = sum(data.Volume[tradeIndeces])
    buyerInitiated = DataFrame(Impact = Vector{Float64}(), NormalizedVolume = Vector{Float64}()); sellerInitiated = DataFrame(Impact = Vector{Float64}(), NormalizedVolume = Vector{Float64}())
    for index in tradeIndeces
        midPriceBeforeTrade = data.MidPrice[index - 1]; midPriceAfterTrade = data.MidPrice[index + 1]
        Δp = log(midPriceAfterTrade) - log(midPriceBeforeTrade)
        ω = (data.Volume[index] / dayVolume) * totalTradeCount
        if !ismissing(Δp) && !ismissing(ω)
            data.Side[index] == -1 ? push!(buyerInitiated, (Δp, ω)) : push!(sellerInitiated, (-Δp, ω))
        end
    end
    return buyerInitiated, sellerInitiated
end
#---------------------------------------------------------------------------------------------------
