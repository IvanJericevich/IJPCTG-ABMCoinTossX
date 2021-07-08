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
    ExtremeLogReturnPercentileDistribution(:TickbyTick, :Upper; lobFile = "Model2L1LOB", format = "png")
    DepthProfile(depthProfile; format = "png")
    InterArrivalTimeDistribution("Model2L1LOB"; format = "png")
    VPIN(data, 50, 10, Millisecond(60 * 30 * 1000))
- TODO: Insert plot annotations for the values of α when fitting power laws and excess kurtosis
- TODO: Change font sizes
- TODO: WHat about mutliple trades in cointyossx l1lob but not in jse
- TODO: SHould comparisons be ontop of each other or seperate
=#
using Distributions, CSV, Plots, DataFrames, StatsPlots, Dates, StatsBase, LaTeXStrings, TimeSeries
clearconsole()
#---------------------------------------------------------------------------------------------------

#----- Generate stylized facts -----#
function StylizedFacts(file1::String, file2::String; resolution = nothing, format::String = "pdf")
    OHLCV("JSEL1LOB", resolution); OHLCV("CoinTossXL1LOB", resolution)
    if isnothing(resolution)
        jsedata = CSV.File(string("Data/JSEL1LOB.csv"), drop = [:MicroPrice, :Spread, :DateTime], missingstring = "missing") |> DataFrame; cointossxdata = CSV.File(string("Data/CoinTossXL1LOB.csv"), drop = [:MicroPrice, :Spread, :DateTime], missingstring = "missing") |> DataFrame
        jselogreturns = diff(log.(filter(x -> !ismissing(x), jsedata[:, :MidPrice]))); cointossxlogreturns = diff(log.(filter(x -> !ismissing(x), cointossxdata[:, :MidPrice])))
    else
        jsedata = CSV.File(string("Data/JSE ", resolution, " Bars.csv"), missingstring = "missing") |> DataFrame; cointossxdata = CSV.File(string("Data/CoinTossX ", resolution, " Bars.csv"), missingstring = "missing") |> DataFrame
        jselogreturns = diff(log.(filter(x -> !ismissing(x), jsedata[:, :MidClose]))); cointossxlogreturns = diff(log.(filter(x -> !ismissing(x), cointossxdata[:, :MidClose])))
    end
    LogReturnDistribution(logreturns; format = format)
    LogReturnAutocorrelation(logreturns, 500; format = format)
    TradeSignAutocorrelation(data, 50; format = format)
    ExtremeLogReturnPercentileDistribution(logreturns; format = format)
    #DepthProfile(depthProfile; format = "png")
end
#---------------------------------------------------------------------------------------------------

#----- Extract OHLCV data -----#
function OHLCV(lobFile::String, resolution)
    l1lob = CSV.File(string("Data/", lobFile, ".csv"), missingstring = "missing", types = Dict(:Type => Symbol)) |> DataFrame
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
                vwap = !isempty(tradesBar) ? sum(tradesBar.TradeVol .* tradesBar.Trade) / sum(tradesBar.TradeVol) : missing
                println(file, string(barTimes[t], ",", midPriceOHLCV, ",", microPriceOHLCV, ",", sum(bar.Volume), ",", vwap))
            end
        end
    end
end
#---------------------------------------------------------------------------------------------------

#----- Log return sample distributions for different time resolutions -----#
function LogReturnDistribution(logreturns::Vector{Float64}; format::String = "pdf")
    normalDistribution = fit(Normal, logreturns)
    abslogreturns = sort(abs.(logreturns)) |> x -> filter(!iszero, x)
    n = length(abslogreturns)
    empiricalcdf = ecdf(logreturns) |> x -> 1 .- x(abslogreturns)
    normalcdf = 1 .- cdf.(normalDistribution, abslogreturns)
    index = minimum([findprev(x -> x > 0, empiricalcdf, n), findprev(x -> x > 0, normalcdf, n)])
    if index != n
        deleteat!(abslogreturns, index:n); deleteat!(empiricalcdf, index:n); deleteat!(normalcdf, index:n)
    end
    distribution = histogram(logreturns, normalize = :pdf, fillcolor = :blue, linecolor = :blue, xlabel = "Log returns", ylabel = "Probability Density", label = "Empirical", legendtitle = "Empirical Distribution", legend = :bottomleft, legendfontsize = 5, legendtitlefontsize = 7, fg_legend = :transparent)
    plot!(distribution, normalDistribution, linecolor = :black, label = "Fitted Normal")
    plot!(distribution, abslogreturns, hcat(empiricalcdf, normalcdf), seriestype = [:scatter :line], markercolor = :blue, markerstrokecolor = :blue, markersize = 3, linecolor = :black, xlabel = "Log returns", ylabel = "Cummulative Density", scale = :log10, legend = false, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, inset = (1, bbox(0.1, 0.1, 0.33, 0.33, :top)), subplot = 2)
    qqplot!(distribution, Normal, logreturns, xlabel = "Normal theoretical quantiles", ylabel = "Sample quantiles", linecolor = :black, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, marker = (:blue, stroke(:blue), 3), legend = false, inset = (1, bbox(0.65, 0.1, 0.33, 0.33, :top)), subplot = 3)
    savefig(distribution, string("Figures/Log-ReturnDistribution.", format))
end
#---------------------------------------------------------------------------------------------------

#----- Log-return and absolute log-return autocorrelation -----#
function LogReturnAutocorrelation(logreturns::Vector{Float64}, lag::Int64; format::String = "pdf")
    autoCorr = autocor(logreturns, 1:lag; demean = false)
    absAutoCorr = autocor(abs.(logreturns), 1:lag; demean = false)
    autoCorrPlot = plot(autoCorr, seriestype = :sticks, linecolor = :black, legend = false, xlabel = "Lag", ylabel = "Autocorrelation")
    plot!(autoCorrPlot, [1.96 / sqrt(length(logreturns)), -1.96 / sqrt(length(logreturns))], seriestype = :hline, line = (:dash, :black, 1))
    plot!(autoCorrPlot, absAutoCorr, seriestype = :sticks, linecolor = :black, legend = false, xlabel = "Lag", ylabel = "Autocorrelation", inset = (1, bbox(0.62, 0.55, 0.33, 0.33, :top)), subplot = 2, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30)
    savefig(autoCorrPlot, "Figures/Log-ReturnAutocorrelation." * format)
end
#---------------------------------------------------------------------------------------------------

#----- Trade sign autocorrealtion -----#
function TradeSignAutocorrelation(data, lag::Int64; format::String = "pdf")
    tradeSigns = filter(x -> x.Type == :Market, data).Side |> y -> map(x -> x == :Buy ? 1 : -1, y)
    autoCorr = autocor(tradeSigns, 1:lag; demean = false)
    autoCorrPlot = plot(autoCorr, seriestype = :sticks, linecolor = :black, legend = false, xlabel = "Lag", ylabel = "Autocorrelation")
    plot!(autoCorrPlot, [quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(tradeSigns)), quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(tradeSigns))], seriestype = :hline, line = (:dash, :black, 1))
    plot!(autoCorrPlot, autoCorr, xscale = :log10, inset = (1, bbox(0.58, 0.0, 0.4, 0.4)), subplot = 2, legend = false, xlabel = "Lag", guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, ylabel = "Autocorrelation", linecolor = :black) #  ", L"(\log_{10})
    savefig(autoCorrPlot, "Figures/Trade-SignAutocorrelation." * format)
end
#---------------------------------------------------------------------------------------------------

#----- Extreme log-return percentile distribution for different time resolutions -----#
function ExtremeLogReturnPercentileDistribution(logreturns::Vector{Float64}; format::String = "pdf")
    upperobservations = logreturns[findall(x -> x >= quantile(logreturns, 0.95), logreturns)]; lowerobservations = -logreturns[findall(x -> x <= quantile(logreturns, 0.05), logreturns)]
    sort!(upperobservations); sort!(lowerobservations)
    upperxₘᵢₙ = minimum(upperobservations); lowerxₘᵢₙ = minimum(lowerobservations)
    upperα = 1 + length(upperobservations) / sum(log.(upperobservations ./ upperxₘᵢₙ)); lowerα = 1 + length(lowerobservations) / sum(log.(lowerobservations ./ lowerxₘᵢₙ))
    upperTheoreticalQuantiles = map(i -> (1 - (i / length(upperobservations))) ^ (-1 / (upperα - 1)) * upperxₘᵢₙ, 1:length(upperobservations)); lowerTheoreticalQuantiles = map(i -> (1 - (i / length(lowerobservations))) ^ (-1 / (lowerα - 1)) * lowerxₘᵢₙ, 1:length(lowerobservations))
    extremePercentileDistributionPlot = plot(upperobservations, seriestype = [:scatter, :line], marker = (:blue, stroke(:blue), :utriangle), normalize = :pdf, linecolor = :blue, xlabel = string("Log return extreme percentiles"), ylabel = "Density", label = ["" "Upper percentiles"], fg_legend = :transparent)#, annotations = (3, y[3], Plots.text(string("α=" α), :left))))
    plot!(extremePercentileDistributionPlot, lowerobservations, seriestype = [:scatter, :line], marker = (:blue, stroke(:blue), :pentagon), normalize = :pdf, linecolor = :blue, label = ["" "Lower percentiles"])
    plot!(extremePercentileDistributionPlot, [upperTheoreticalQuantiles upperTheoreticalQuantiles], [upperobservations upperTheoreticalQuantiles], seriestype = [:scatter :line], inset = (1, bbox(0.2, 0.03, 0.34, 0.34, :top)), subplot = 2, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, legend = :none, xlabel = "Power-Law Theoretical Quantiles", ylabel = "Sample Quantiles", linecolor = :black, markercolor = :blue, markerstrokecolor = :blue, markershape = :utriangle, markersize = 3, fg_legend = :transparent)
    plot!(extremePercentileDistributionPlot, [lowerTheoreticalQuantiles lowerTheoreticalQuantiles], [lowerobservations lowerTheoreticalQuantiles], seriestype = [:scatter :line], subplot = 2, linecolor = :black, markercolor = :blue, markerstrokecolor = :blue, markershape = :pentagon, markersize = 3)
    savefig(extremePercentileDistributionPlot, string("Figures/ExtremeLog-ReturnPercentilesDistribution.", format))
end
#---------------------------------------------------------------------------------------------------

#----- Depth profile -----#
function DepthProfile(profile::Array{Union{Missing, Int64}, 2}; format::String = "pdf")
    μ = map(i -> mean(skipmissing(profile[:, i])), 1:size(profile, 2))
    depthProfile = plot(-(1:7), μ[1:7], seriestype = [:scatter, :line], marker = (:blue, stroke(:blue), :utriangle), linecolor = :blue, label = ["" "Bid profile"], xlabel = "Price level of limit orders (<0: bids; >0: asks)", ylabel = "Volume", fg_legend = :transparent)
    plot!(depthProfile, 1:7, μ[8:14], seriestype = [:scatter, :line], marker = (:red, stroke(:red), :dtriangle), linecolor = :red, label = ["" "Ask profile"])
    savefig(depthProfile, string("Figures/DepthProfile.", format))
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
