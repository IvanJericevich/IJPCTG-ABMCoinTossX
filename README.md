# Simulation and estimation of an agent-based market-model with a matching engine

## Authors
* Ivan Jericevich
* Patrick Chang
* Tim Gebbie


## Description
We present an agent-based model in an artificial exchange framework that is able to recover price-impact while replicating the conventional stylised facts by directly interacting with a decoupled Matching Engine (ME) that manages orders away from the software relating to agent interactions and rules. This repo contains scripts for the replication of results in the pubished article. The scripts in this repository provide a number of different functions:
* A flexible asynchronous ABM for interacting with CoinTossX
* A Julia API for CoinTossX along with the functionality for polling the market-data feed
* Code for the recovery of sensitivity analysis and calibration results
* Methodology for the replication of many stylized facts (including price impact)
* A simple implementation of the Nelder-Mead with Threshold Accepting optimization heuristic
* Construction of the method of moments objective and implementation of the simulated minimum distance method in the context of financial market ABMs
* JSE and CoinTossX data cleaning functions

The key innovation here is that the model updates no longer occur in chronological or calendar time. Here the Model Management System (MMS) polls the market data feed based on asynchronous model time, rather than be entirely driven by the event time of the market. This framework provides a basis for future development.

## Prerequisites
* [Julia](https://julialang.org) programming langauge
* A text editor or IDE such as [Atom](https://flight-manual.atom.io/getting-started/sections/installing-atom/), [VS Code](https://code.visualstudio.com/download) or [Jupyter](https://jupyter.org/install)
* CoinTossX matching engine. The version of CoinTossX used is a fork from the official release of version [v1.1.0](https://github.com/dharmeshsing/CoinTossX/tree/v1.1.0) and contains added functionalty necessary for ABM simulation and calibration. This modified version can be found [here](https://github.com/IvanJericevich/CoinTossX/tree/abm-experiment). For instructions on how to install, start, and use CoinTossX refer to [CoinTossX](https://github.com/dharmeshsing/CoinTossX)

## Usage
Clone the repository
```sh
git clone https://github.com/IvanJericevich/IJPCTG-ABMCoinTossX.git
```
Packages can be installed using
```julia
Pkg.add("...")
```
Each script contains examples of how to run the relevant functions. [CoinTossXUtilities.jl](Scripts/CoinTossXUtilities.jl) contains the important functions for interacting with CoinTossX. [ABMVolatilityAuctionProxy.jl](Scripts/ABMVolatilityAuctionProxy.jl) is the script used to simulate using the ABM and submit generated orders via the model to CoinTossX. [Calibration.jl](Scripts/Calibration.jl) calibrates parameters to empirical JSE TAQ data. [SensitivityAnalysis.jl](Scripts/SensitivityAnalysis.jl) provides code for the recovery of sensitivity analysis results and figures. Important to the calibration of parameters and sensitivity analysis is the algorithm for cleaning and re-formatting raw data from both COinTossX and JSE. This is demonstrated in [DataCleaning](DataCleaning).

### Data structure
Provided in the [Data](Data) folder is cleaned level 1 LOB data from calibrated parameters applied to CoinTossX (for visualisation), along with the calibrated parameters (with confidence intervals). Additionally, the NMTA optimization stacktrace with the bootstrapped covariance matrix of empirical log-returns are provided as well. The exact files and their contents are as follows:
* [Data/Trader.csv](Data/Trades_1.csv) contains the trader menmonics of all trader agent types.
* [Data/ParametersCI.csv](Data/ParametersCI.csv) contains the error bars around calibrated parameters.
* [Data/MomentsCI.csv](Data/MomentsCI.csv) contains the error bars around moments computed from simulationds with calibrated parameters.
* [Data/Calibration/L1LOB.csv](Data/Calibration/L1LOB.csv) contains the cleaned level 1 LOB data from simulations with calibrated parameters.
* [Data/Calibration/OptimisationResults.jld](Data/Calibration/OptimisationResults.jld) is the data file for the NMTA optimization stacktrace.
* [Data/Calibration/W.jld](Data/Calibration/W.jld) is the data file for the bootstrapped covariance matrix of empirical log-returns.
