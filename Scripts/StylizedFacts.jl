#=
StylizedFacts:
- Julia version: 1.5.3
- Authors: Ivan Jericevich, Patrick Chang, Tim Gebbie
- Function: Plot the stylized facts of HFT data for different time resolutions
- Structure:
    1. Generate stylized facts
    2. Extract OHLCV data
    3. Log return sample distributions for different time resolutions
    4. Log-return and absolute log-return autocorrelation
    5. Trade sign autocorrealtion
    6. Extreme log-return percentile distribution for different time resolutions
    7. Volume-volatility correlation
    8. Depth profile
    9. Price Impact
- Examples
    PriceImpact("Sensitivity"); PriceImpact("JSE")
    StylizedFacts("CoinTossX"); StylizedFacts("JSE")
=#
using Distributions, CSV, Plots, StatsPlots, Dates, StatsBase, Distributed, DataFrames, Distributed
addprocs(2)
@everywhere import Statistics.var
@everywhere using DataFrames
clearconsole()
#---------------------------------------------------------------------------------------------------

#----- Generate stylized facts -----#
function StylizedFacts(exchange::String; format::String = "pdf")
    data = CSV.File(string("Data/", exchange, "/L1LOB.csv"), drop = [:MidPrice, :Spread, :Price, :Volume], missingstring = "missing", types = Dict(:Type => Symbol)) |> DataFrame
    if exchange == "JSE"
        filter!(y -> y.DateTime < DateTime("2019-01-04T09:00:41"), data)
    end
    data.Date = Date.(data.DateTime)
    uniqueDays = unique(data.Date)
    logreturns = map(day -> diff(log.(skipmissing(data[searchsorted(data.Date, day), :MicroPrice]))), uniqueDays) |> x -> reduce(vcat, x)
    LogReturnDistribution(exchange, logreturns; format = format)
    LogReturnAutocorrelation(exchange, logreturns, 500; format = format)
    TradeSignAutocorrelation(exchange, data, 500; format = format)
    ExtremeLogReturnPercentileDistribution(exchange, logreturns; format = format)
    VolumeVolatilityCorrelation(exchange, data; N = 500, format = format)
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
function LogReturnDistribution(exchange::String, logreturns::Vector{Float64}; format::String = "pdf")
    color = exchange == "CoinTossX" ? :green : (exchange == "JSE" ? :purple : :orange)
    NormalDistribution = fit(Normal, logreturns)
    distribution = histogram(logreturns, normalize = :pdf, fillcolor = color, linecolor = color, xlabel = "Log returns", ylabel = "Probability Density", label = "Empirical", legendtitle = "Distribution", legend = :topright, legendfontsize = 5, legendtitlefontsize = 7, fg_legend = :transparent, ylim = (0, 3000), xlim = (-0.012, 0.01))
    plot!(distribution, NormalDistribution, line = (:black, 2), label = "Fitted Normal")
    qqplot!(distribution, Normal, logreturns, xlabel = "Normal theoretical quantiles", ylabel = "Sample quantiles", linecolor = :black, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, marker = (color, stroke(color), 3), legend = false, inset = (1,bbox(0.15, 0.03, 0.4, 0.4)), subplot = 2, title = "Normal QQ-plot", titlefontsize = 7)
    savefig(distribution, string("Figures/", exchange, "Log-ReturnDistribution.", format))
end
#---------------------------------------------------------------------------------------------------

#----- Log-return and absolute log-return autocorrelation -----#
function LogReturnAutocorrelation(exchange::String, logreturns::Vector{Float64}, lag::Int64; format::String = "pdf")
    color = exchange == "CoinTossX" ? :green : (exchange == "JSE" ? :purple : :orange)
    autoCorr = autocor(logreturns, 1:lag; demean = false)
    absAutoCorr = autocor(abs.(logreturns), 1:lag; demean = false)
    autoCorrPlot = plot(autoCorr, seriestype = [:sticks, :scatter], marker = (color, stroke(color), 3), linecolor = :black, xlabel = "Lag", ylabel = "Autocorrelation", legend = false, ylim = (-0.4, 0.2))
    plot!(autoCorrPlot, [1.96 / sqrt(length(logreturns)), -1.96 / sqrt(length(logreturns))], seriestype = :hline, line = (:dash, :black, 2))
    plot!(autoCorrPlot, absAutoCorr, seriestype = :scatter, marker = (color, stroke(color), 3), legend = false, xlabel = "Lag", ylabel = "Autocorrelation", inset = (1, bbox(0.62, 0.5, 0.4, 0.4, :top)), subplot = 2, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, title = "Absolute log-return autocorrelation", titlefontsize = 7)
    savefig(autoCorrPlot, string("Figures/", exchange, "Log-ReturnAutocorrelation.", format))
end
#---------------------------------------------------------------------------------------------------

#----- Trade sign autocorrealtion -----#
function TradeSignAutocorrelation(exchange::String, data::DataFrame, lag::Int64; format::String = "pdf")
    color = exchange == "CoinTossX" ? :green : (exchange == "JSE" ? :purple : :orange)
    tradeSigns = data[findall(x -> x == :Market, data.Type), :Side]
    autoCorr = autocor(tradeSigns, 1:lag; demean = false)
    autoCorrPlot = plot(autoCorr, seriestype = :scatter, linecolor = :black, marker = (color, stroke(color), 3), legend = false, xlabel = "Lag", ylabel = "Autocorrelation", fg_legend = :transparent, ylim = (-0.1, 0.8))
    plot!(autoCorrPlot, [quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(tradeSigns)), quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(tradeSigns))], seriestype = :hline, line = (:dash, :black, 2))
    plot!(autoCorrPlot, autoCorr, xscale = :log10, inset = (1, bbox(0.58, 0.1, 0.4, 0.4)), subplot = 2, legend = false, xlabel = "Lag", guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, ylabel = "Autocorrelation", linecolor = color, title = "Log-scale order-flow autocorrelation", titlefontsize = 7) #  ", L"(\log_{10})
    savefig(autoCorrPlot, string("Figures/", exchange, "Trade-SignAutocorrelation.", format))
end
#---------------------------------------------------------------------------------------------------

#----- Extreme log-return percentile distribution for different time resolutions -----#
function ExtremeLogReturnPercentileDistribution(exchange::String, logreturns::Vector{Float64}; format::String = "pdf")
    color = exchange == "CoinTossX" ? :green : (exchange == "JSE" ? :purple : :orange)
    upperobservations = logreturns[findall(x -> x >= quantile(logreturns, 0.95), logreturns)]; lowerobservations = -logreturns[findall(x -> x <= quantile(logreturns, 0.05), logreturns)]
    sort!(upperobservations); sort!(lowerobservations)
    upperxₘᵢₙ = minimum(upperobservations); lowerxₘᵢₙ = minimum(lowerobservations)
    upperα = 1 + length(upperobservations) / sum(log.(upperobservations ./ upperxₘᵢₙ)); lowerα = 1 + length(lowerobservations) / sum(log.(lowerobservations ./ lowerxₘᵢₙ))
    upperTheoreticalQuantiles = map(i -> (1 - (i / length(upperobservations))) ^ (-1 / (upperα - 1)) * upperxₘᵢₙ, 1:length(upperobservations)); lowerTheoreticalQuantiles = map(i -> (1 - (i / length(lowerobservations))) ^ (-1 / (lowerα - 1)) * lowerxₘᵢₙ, 1:length(lowerobservations))
    extremePercentileDistributionPlot = density(upperobservations, seriestype = [:scatter, :line], marker = (color, stroke(color), :utriangle), linecolor = color, xlabel = string("Log return extreme percentiles"), ylabel = "Density", label = string("Upper percentiles - α = ", round(upperα, digits = 3)), legend = :topright, fg_legend = :transparent)
    density!(extremePercentileDistributionPlot, lowerobservations, seriestype = [:scatter, :line], marker = (color, stroke(color), :dtriangle), linecolor = color, label = string("Lower percentiles - α = ", round(lowerα, digits = 3)))
    plot!(extremePercentileDistributionPlot, hcat(upperTheoreticalQuantiles, upperTheoreticalQuantiles), hcat(upperobservations, upperTheoreticalQuantiles), scale = :log10, seriestype = [:scatter :line], inset = (1, bbox(0.2, 0.03, 0.34, 0.34, :top)), subplot = 2, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, legend = :none, xlabel = "Power-Law Theoretical Quantiles", ylabel = "Sample Quantiles", linecolor = :black, marker = (color, stroke(color), 3, [:utriangle :none]), fg_legend = :transparent, title = "Power-Law QQ-plot", titlefontsize = 7)
    plot!(extremePercentileDistributionPlot, [lowerTheoreticalQuantiles lowerTheoreticalQuantiles], [lowerobservations lowerTheoreticalQuantiles], seriestype = [:scatter :line], subplot = 2, linecolor = :black, marker = (color, stroke(color), 3, [:dtriangle :none]))
    savefig(extremePercentileDistributionPlot, string("Figures/", exchange, "ExtremeLog-ReturnPercentilesDistribution.", format))
end
#---------------------------------------------------------------------------------------------------

#----- Volume-volatility correlation -----#
function VolumeVolatilityCorrelation(exchange::String, data::DataFrame; N = 5000, format::String = "pdf")
    tradeIndeces = findall(x -> x == :Market, data.Type)
    days = unique(data.Date)
    variances = @distributed (hcat) for day in days
        dayIndeces = tradeIndeces[searchsorted(data[tradeIndeces, :Date], day)]
        σ = map(i -> var(diff(log.(skipmissing(data[1:(dayIndeces[i]), :MicroPrice])))), 1:N)
        σ
    end
    color = exchange == "CoinTossX" ? :green : :purple
    correlation = plot(1:N, mean(variances, dims = 2), seriestype = :line, linecolor = color, xlabel = "Variance", ylabel = "Number of trades", legend = false, scale = :log10)
    savefig(correlation, string("Figures/", exchange, "Volume-VolatilityCorrelation.", format))
end
#---------------------------------------------------------------------------------------------------

#----- Depth profile -----#
function DepthProfile(profile::Array{Union{Missing, Int64}, 2}; format::String = "pdf")
    μ = map(i -> mean(skipmissing(profile[:, i])), 1:size(profile, 2))
    depthProfile = plot(-(1:7), μ[1:7], seriestype = [:scatter, :line], marker = (:blue, stroke(:blue), :utriangle), linecolor = :blue, label = ["" "Bid profile"], xlabel = "Price level of limit orders (<0: bids; >0: asks)", ylabel = "Volume", fg_legend = :transparent)#, yscale = :log10
    plot!(depthProfile, 1:7, μ[8:14], seriestype = [:scatter, :line], marker = (:red, stroke(:red), :dtriangle), linecolor = :red, label = ["" "Ask profile"])
    savefig(depthProfile, string("Figures/DepthProfile.", format))
end
#---------------------------------------------------------------------------------------------------

#----- Extract price-impact data -----#
function PriceImpact(exchange::String; format::String = "pdf")
    data = CSV.File(string("Data/", exchange, "/L1LOB.csv"), drop = [:MicroPrice, :Price, :Spread], missingstring = "missing", types = Dict(:Type => Symbol)) |> DataFrame
    buyerInitiated = DataFrame(Impact = Vector{Float64}(), NormalizedVolume = Vector{Float64}()); sellerInitiated = DataFrame(Impact = Vector{Float64}(), NormalizedVolume = Vector{Float64}())
    days = unique(x -> Date(x), data.DateTime)
    tradeIndeces = findall(x -> x == :Market, data.Type)
    totalTradeCount = length(tradeIndeces)
    for day in days
        dayTradeIndeces = tradeIndeces[searchsorted(data[tradeIndeces, :DateTime], day, by = Date)]
        dayVolume = sum(data[dayTradeIndeces, :Volume])
        for index in dayTradeIndeces
            midPriceBeforeTrade = data[index - 1, :MidPrice]; midPriceAfterTrade = data[index + 1, :MidPrice]
            Δp = log(midPriceAfterTrade) - log(midPriceBeforeTrade)
            ω = (data[index, :Volume] / dayVolume) * (totalTradeCount / length(days))
            if !ismissing(Δp) && !ismissing(ω)
                data[index, :Side] == 1 ? push!(buyerInitiated, (Δp, ω)) : push!(sellerInitiated, (-Δp, ω))
            end
        end
    end
    filter!(x -> !isnan(x.Impact) && !isnan(x.NormalizedVolume), buyerInitiated); filter!(x -> !isnan(x.Impact) && !isnan(x.NormalizedVolume), sellerInitiated)
    normalisedVolumeBins = 10 .^ (range(-1, 1, length = 21))
    Δp = fill(NaN, (length(normalisedVolumeBins), 2)); ω = fill(NaN, (length(normalisedVolumeBins), 2)) # Column 1 is buy; column 2 is sell
    for i in 2:length(normalisedVolumeBins)
        binIndeces = (findall(x -> normalisedVolumeBins[i - 1] < x <= normalisedVolumeBins[i], buyerInitiated.NormalizedVolume), findall(x -> normalisedVolumeBins[i - 1] < x <= normalisedVolumeBins[i], sellerInitiated.NormalizedVolume))
        if !isempty(first(binIndeces))
            Δp[i - 1, 1] = mean(buyerInitiated[first(binIndeces), :Impact]); ω[i, 1] = mean(buyerInitiated[first(binIndeces), :NormalizedVolume])
        end
        if !isempty(last(binIndeces))
            Δp[i - 1, 2] = mean(sellerInitiated[last(binIndeces), :Impact]); ω[i, 2] = mean(sellerInitiated[last(binIndeces), :NormalizedVolume])
        end
    end
    indeces = findall(vec(any(x -> !isnan(x), Δp, dims = 2) .* any(x -> !isnan(x), ω, dims = 2)))
    priceImpact = plot(ω[2:(end-3), :], Δp[2:(end-3), :], scale = :log10, seriestype = [:scatter, :line], markershape = [:utriangle :dtriangle], markercolor = [:blue :red], markerstrokecolor = [:blue :red], markersize = 3, linecolor = [:blue :red], xlabel = "ω*", ylabel = "Δp*", label = ["" "" "Buyer initiated" "Seller initiated"], legend = :topleft, fg_legend = :transparent)
    savefig(priceImpact, string("Figures/", exchange, "PriceImpact.", format))
end
#---------------------------------------------------------------------------------------------------
