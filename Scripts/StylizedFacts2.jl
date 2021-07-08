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
function StylizedFacts(; resolution = nothing, format::String = "pdf")
    if isnothing(resolution)
        jsedata = CSV.File(string("Data/JSEL1LOB.csv"), drop = [:MidPrice, :Spread, :DateTime], missingstring = "missing") |> DataFrame; cointossxdata = CSV.File(string("Data/CoinTossXL1LOB.csv"), drop = [:MidPrice, :Spread, :DateTime], missingstring = "missing") |> DataFrame
        jselogreturns = diff(log.(filter(x -> !ismissing(x), jsedata[:, :MicroPrice]))); cointossxlogreturns = diff(log.(filter(x -> !ismissing(x), cointossxdata[:, :MicroPrice])))
    else
        OHLCV("JSEL1LOB", resolution); OHLCV("CoinTossXL1LOB", resolution)
        jsedata = CSV.File(string("Data/JSE ", resolution, " Bars.csv"), missingstring = "missing") |> DataFrame; cointossxdata = CSV.File(string("Data/CoinTossX ", resolution, " Bars.csv"), missingstring = "missing") |> DataFrame
        jselogreturns = diff(log.(filter(x -> !ismissing(x), jsedata[:, :MidClose]))); cointossxlogreturns = diff(log.(filter(x -> !ismissing(x), cointossxdata[:, :MidClose])))
    end
    LogReturnDistribution(jselogreturns, cointossxlogreturns; format = format)
    LogReturnAutocorrelation(jselogreturns, cointossxlogreturns, 500; format = format)
    TradeSignAutocorrelation(jsedata, cointossxdata, 50; format = format)
    ExtremeLogReturnPercentileDistribution(jselogreturns, cointossxlogreturns; format = format)
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
jsedata = CSV.File(string("Data/JSEL1LOB.csv"), drop = [:MidPrice, :Spread, :DateTime], types = Dict(:Type => Symbol), missingstring = "missing") |> DataFrame
jselogreturns1 = diff(log.(filter(x -> !ismissing(x), jsedata[1:81309, :MicroPrice])))
jselogreturns2 = diff(log.(filter(x -> !ismissing(x), jsedata[81310:end, :MicroPrice])))


ExtremeLogReturnPercentileDistribution(jselogreturns1, jselogreturns2; format = "png")
#----- Log return sample distributions for different time resolutions -----#
function LogReturnDistribution(jselogreturns::Vector{Float64}, cointossxlogreturns::Vector{Float64}; format::String = "pdf")
    jseNormalDistribution = fit(Normal, jselogreturns); cointossxNormalDistribution = fit(Normal, cointossxlogreturns)
    jseabslogreturns = sort(abs.(jselogreturns)) |> x -> filter(!iszero, x); cointossxabslogreturns = sort(abs.(cointossxlogreturns)) |> x -> filter(!iszero, x)
    jsen = length(jseabslogreturns); cointossxn = length(cointossxabslogreturns)
    jseempiricalcdf = ecdf(jselogreturns) |> x -> 1 .- x(jseabslogreturns); cointossxempiricalcdf = ecdf(cointossxlogreturns) |> x -> 1 .- x(cointossxabslogreturns)
    jsenormalcdf = 1 .- cdf.(jseNormalDistribution, jseabslogreturns); cointossxnormalcdf = 1 .- cdf.(cointossxNormalDistribution, cointossxabslogreturns)
    jseindex = minimum([findprev(x -> x > 0, jseempiricalcdf, jsen), findprev(x -> x > 0, jsenormalcdf, jsen)]); cointossxindex = minimum([findprev(x -> x > 0, cointossxempiricalcdf, cointossxn), findprev(x -> x > 0, cointossxnormalcdf, cointossxn)])
    if jseindex != jsen deleteat!(jseabslogreturns, jseindex:jsen); deleteat!(jseempiricalcdf, jseindex:jsen); deleteat!(jsenormalcdf, jseindex:jsen) end
    if cointossxindex != cointossxn deleteat!(cointossxabslogreturns, cointossxindex:cointossxn); deleteat!(cointossxempiricalcdf, cointossxindex:cointossxn); deleteat!(cointossxnormalcdf, cointossxindex:cointossxn) end
    distribution = histogram(jselogreturns, normalize = :pdf, fill = (:purple, 0.2), line = (:purple, 0.2), xlabel = "Log returns", ylabel = "Probability Density", label = "JSE - Empirical", legendtitle = "Distribution", legend = :bottomleft, legendfontsize = 5, legendtitlefontsize = 7, fg_legend = :transparent)
    histogram!(distribution, cointossxlogreturns, normalize = :pdf, fill = (:green, 0.2), line = (:green, 0.2), label = "CoinTossX - Empirical")
    plot!(distribution, jseNormalDistribution, line = (:purple, 2), label = "JSE - Fitted Normal")
    plot!(distribution, cointossxNormalDistribution, line = (:green, 2), label = "CoinTossX - Fitted Normal")
    plot!(distribution, jseabslogreturns, hcat(jseempiricalcdf, jsenormalcdf), seriestype = [:scatter :line], markercolor = :purple, markerstrokecolor = :purple, markersize = 3, linecolor = :black, xlabel = "Log returns", ylabel = "Cummulative Density", scale = :log10, legend = false, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, inset = (1, bbox(0.1, 0.1, 0.33, 0.33, :top)), subplot = 2)
    plot!(distribution, cointossxabslogreturns, hcat(cointossxempiricalcdf, cointossxnormalcdf), seriestype = [:scatter :line], markercolor = :green, markerstrokecolor = :green, markersize = 3, linecolor = :black, scale = :log10, subplot = 2)
    qqplot!(distribution, Normal, jselogreturns, xlabel = "Normal theoretical quantiles", ylabel = "Sample quantiles", linecolor = :black, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, marker = (:purple, stroke(:purple), 3), legend = false, inset = (1, bbox(0.65, 0.1, 0.33, 0.33, :top)), subplot = 3)
    qqplot!(distribution, Normal, cointossxlogreturns, linecolor = :black, marker = (:green, stroke(:green), 3), subplot = 3)
    savefig(distribution, string("Figures/Log-ReturnDistribution.", format))
end
#---------------------------------------------------------------------------------------------------

#----- Log-return and absolute log-return autocorrelation -----#
function LogReturnAutocorrelation(jselogreturns::Vector{Float64}, cointossxlogreturns::Vector{Float64}, lag::Int64; format::String = "pdf")
    autoCorr = hcat(autocor(jselogreturns, 1:lag; demean = false), autocor(cointossxlogreturns, 1:lag; demean = false))
    absAutoCorr = hcat(autocor(abs.(jselogreturns), 1:lag; demean = false), autocor(abs.(cointossxlogreturns), 1:lag; demean = false))
    autoCorrPlot = plot(autoCorr, seriestype = [:sticks, :scatter], markercolor = [:purple :green], markerstrokecolor = [:purple :green], markersize = 3, linecolor = :black, label = ["" "" "JSE" "CoinTossX"], xlabel = "Lag", ylabel = "Autocorrelation", legend = :bottomleft)
    plot!(autoCorrPlot, hcat([1.96 / sqrt(length(jselogreturns)), -1.96 / sqrt(length(jselogreturns))], [1.96 / sqrt(length(cointossxlogreturns)), -1.96 / sqrt(length(cointossxlogreturns))]), seriestype = :hline, line = (:dash, [:purple :green], 2), label = "")
    plot!(autoCorrPlot, absAutoCorr, seriestype = [:sticks, :scatter], linecolor = :black, markercolor = [:purple :green], markerstrokecolor = [:purple :green], markersize = 3, legend = false, xlabel = "Lag", ylabel = "Autocorrelation", inset = (1, bbox(0.62, 0.55, 0.33, 0.33, :top)), subplot = 2, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30)
    savefig(autoCorrPlot, "Figures/Log-ReturnAutocorrelation." * format)
end
#---------------------------------------------------------------------------------------------------

#----- Trade sign autocorrealtion -----#
function TradeSignAutocorrelation(jsedata::DataFrame, cointossxdata::DataFrame, lag::Int64; format::String = "pdf")
    function TradeSigns(data::DataFrame)
        tradeSigns = Vector{Int64}(); index = 0
        for i in 1:nrow(data)
            if i > index && data[i, :Type] == :Market
                push!(tradeSigns, data[i, :Side])
                index = findnext(x -> x != :Market, data.Type, i)
            end
        end
        return tradeSigns
    end
    jseTradeSigns = TradeSigns(jsedata); cointossxTradeSigns = TradeSigns(cointossxdata)
    autoCorr = hcat(autocor(jseTradeSigns, 1:lag; demean = false), autocor(cointossxTradeSigns, 1:lag; demean = false))
    autoCorrPlot = plot(autoCorr, seriestype = [:sticks, :scatter], linecolor = :black, marker = ([:purple :green], 3), markerstrokecolor = [:purple :green], label = ["" "" "JSE" "CoinTossX"], legend = :topleft, xlabel = "Lag", ylabel = "Autocorrelation")
    plot!(autoCorrPlot, hcat([quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(jseTradeSigns)), quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(jseTradeSigns))], [quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(cointossxTradeSigns)), quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(cointossxTradeSigns))]), seriestype = :hline, line = (:dash, [:purple :green], 2), label = "")
    plot!(autoCorrPlot, autoCorr, xscale = :log10, inset = (1, bbox(0.58, 0.0, 0.4, 0.4)), subplot = 2, legend = false, xlabel = "Lag", guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, ylabel = "Autocorrelation", linecolor = [:purple :green]) #  ", L"(\log_{10})
    savefig(autoCorrPlot, "Figures/Trade-SignAutocorrelation." * format)
end
#---------------------------------------------------------------------------------------------------

#----- Extreme log-return percentile distribution for different time resolutions -----#
function ExtremeLogReturnPercentileDistribution(jselogreturns::Vector{Float64}, cointossxlogreturns::Vector{Float64}; format::String = "pdf")
    jseupperobservations = jselogreturns[findall(x -> x >= quantile(jselogreturns, 0.95), jselogreturns)]; jselowerobservations = -jselogreturns[findall(x -> x <= quantile(jselogreturns, 0.05), jselogreturns)]
    cointossxupperobservations = cointossxlogreturns[findall(x -> x >= quantile(cointossxlogreturns, 0.95), cointossxlogreturns)]; cointossxlowerobservations = -cointossxlogreturns[findall(x -> x <= quantile(cointossxlogreturns, 0.05), cointossxlogreturns)]
    sort!(jseupperobservations); sort!(jselowerobservations); sort!(cointossxupperobservations); sort!(cointossxlowerobservations)
    jseupperxₘᵢₙ = minimum(jseupperobservations); jselowerxₘᵢₙ = minimum(jselowerobservations); cointossxupperxₘᵢₙ = minimum(cointossxupperobservations); cointossxlowerxₘᵢₙ = minimum(cointossxlowerobservations)
    jseupperα = 1 + length(jseupperobservations) / sum(log.(jseupperobservations ./ jseupperxₘᵢₙ)); jselowerα = 1 + length(jselowerobservations) / sum(log.(jselowerobservations ./ jselowerxₘᵢₙ))
    cointossxupperα = 1 + length(cointossxupperobservations) / sum(log.(cointossxupperobservations ./ cointossxupperxₘᵢₙ)); cointossxlowerα = 1 + length(cointossxlowerobservations) / sum(log.(cointossxlowerobservations ./ cointossxlowerxₘᵢₙ))
    jseUpperTheoreticalQuantiles = map(i -> (1 - (i / length(jseupperobservations))) ^ (-1 / (jseupperα - 1)) * jseupperxₘᵢₙ, 1:length(jseupperobservations)); jseLowerTheoreticalQuantiles = map(i -> (1 - (i / length(jselowerobservations))) ^ (-1 / (jselowerα - 1)) * jselowerxₘᵢₙ, 1:length(jselowerobservations))
    cointossxUpperTheoreticalQuantiles = map(i -> (1 - (i / length(cointossxupperobservations))) ^ (-1 / (cointossxupperα - 1)) * cointossxupperxₘᵢₙ, 1:length(cointossxupperobservations)); cointossxLowerTheoreticalQuantiles = map(i -> (1 - (i / length(cointossxlowerobservations))) ^ (-1 / (cointossxlowerα - 1)) * cointossxlowerxₘᵢₙ, 1:length(cointossxlowerobservations))
    extremePercentileDistributionPlot = plot(jseupperobservations, seriestype = [:scatter, :line], marker = (:purple, stroke(:purple), :utriangle), normalize = :pdf, linecolor = :purple, xlabel = string("Log return extreme percentiles"), ylabel = "Density", label = ["" string("JSE upper percentiles - α = ", round(jseupperα, digits = 3))], legend = :bottomleft, fg_legend = :transparent)
    plot!(extremePercentileDistributionPlot, jselowerobservations, seriestype = [:scatter, :line], marker = (:purple, stroke(:purple), :dtriangle), normalize = :pdf, linecolor = :purple, label = ["" string("JSE lower percentiles - α = ", round(jselowerα, digits = 3))])
    plot!(extremePercentileDistributionPlot, cointossxupperobservations, seriestype = [:scatter, :line], marker = (:green, stroke(:green), :utriangle), normalize = :pdf, linecolor = :green, label = ["" string("CoinTossX upper percentiles - α = ", round(cointossxupperα, digits = 3))])
    plot!(extremePercentileDistributionPlot, cointossxlowerobservations, seriestype = [:scatter, :line], marker = (:green, stroke(:green), :dtriangle), normalize = :pdf, linecolor = :green, label = ["" string("CoinTossX lower percentiles - α = ", round(cointossxlowerα, digits = 3))])
    plot!(extremePercentileDistributionPlot, hcat(jseUpperTheoreticalQuantiles, jseUpperTheoreticalQuantiles), hcat(jseupperobservations, jseUpperTheoreticalQuantiles), seriestype = [:scatter :line], inset = (1, bbox(0.2, 0.03, 0.34, 0.34, :top)), subplot = 2, guidefontsize = 7, tickfontsize = 5, xrotation = 30, yrotation = 30, legend = :none, xlabel = "Power-Law Theoretical Quantiles", ylabel = "Sample Quantiles", linecolor = :black, markercolor = :purple, markerstrokecolor = :purple, markershape = :utriangle, markersize = 3, fg_legend = :transparent)
    plot!(extremePercentileDistributionPlot, hcat(jseLowerTheoreticalQuantiles, jseLowerTheoreticalQuantiles), hcat(jselowerobservations, jseLowerTheoreticalQuantiles), seriestype = [:scatter :line], subplot = 2, linecolor = :black, markercolor = :purple, markerstrokecolor = :purple, markershape = :dtriangle, markersize = 3)
    plot!(extremePercentileDistributionPlot, hcat(cointossxUpperTheoreticalQuantiles, cointossxUpperTheoreticalQuantiles), hcat(cointossxupperobservations, cointossxUpperTheoreticalQuantiles), seriestype = [:scatter :line], subplot = 2, linecolor = :black, markercolor = :green, markerstrokecolor = :green, markershape = :utriangle, markersize = 3, fg_legend = :transparent)
    plot!(extremePercentileDistributionPlot, hcat(cointossxLowerTheoreticalQuantiles, cointossxLowerTheoreticalQuantiles), hcat(cointossxlowerobservations, cointossxLowerTheoreticalQuantiles), seriestype = [:scatter :line], subplot = 2, linecolor = :black, markercolor = :green, markerstrokecolor = :green, markershape = :dtriangle, markersize = 3)
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
